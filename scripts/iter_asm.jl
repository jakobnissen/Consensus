"""
    SelectTemplates

Module containing code to select templates after the initial round of mapping
"""
module SelectTemplates

using Influenza: Segment, Reference, load_references
using FASTX: FASTA
using BioSequences: LongDNASeq
import KMATools: KMATools

imap(f) = x -> Iterators.map(f, x)

# This is not optimized for speed, but I seriously doubt it's going to be a problem
"""
    read_res(path, segment_map) => Vector{String}

Parse the kma.res file at `path`, returning a vector of names of the templates
that should be used.
"""
function get_templates(
    res_path::AbstractString,
    ref_dict::Dict{String, Reference} # name => ref
)::Vector{String}
    rows = open(io -> KMATools.parse_res(io, res_path), res_path)
    bysegment = Dict{Segment, Vector{eltype(rows)}}()
    for row in rows
        push!(get!(valtype(bysegment), bysegment, ref_dict[row.template].segment), row)
    end

    # Criterion by which to select the best segment
    function isworse(x, than)
        # Coverage is most critical: Without it, it's impossible to even make an
        # alignment in the first place. But we can't expect more than around 90%
        # due to deletions, or spurious termini in references etc.
        if min(x.tcov, than.tcov) < 0.9
            return x.tcov < than.tcov
        # Why is identity important? Not sure it needs to be here, but a bad ref
        # could get higher depth due to crossmapping from another segment, e.g. at
        # termini. This factor is small if identity is sufficiently high
        elseif min(x.tid, than.tid) < 0.85
            return x.tid < than.tid
        # If coverage and id is OK, depth is the most important factor, since
        # this is what determines what the majority (consensus) variant is.
        else
            return x.depth < than.depth
        end
    end

    best = Dict(k => partialsort!(v, 1, lt=isworse, rev=true) for (k, v) in bysegment)
    # For the others (non-best), we have different criteria: We pick all hits
    # with id > 85%, cov > 90% and depth > 0.025 * highest_depth
    highest_depth = Dict(k => first(sort(v, by=i -> i.depth, rev=true)) for (k,v) in bysegment)

    result = [b.template for b in values(best)]
    for (segment, v) in bysegment
        vfilt = filter(v[2:end]) do row # exclude best segment at index 1
            row.tcov ≥ 0.90 &&
            row.tid ≥ 0.85 &&
            row.depth > max(10.0, 0.025 * highest_depth[segment].depth)
        end
        append!(result, [i.template for i in vfilt])
    end
    return result
end

function dump_sequences(
    out_path::AbstractString,
    headers::Vector{String},
    ref_dict::Dict{String, Reference}
)::Dict{String, LongDNASeq} # map of name => seq of all templates
    open(FASTA.Writer, out_path) do writer
        result = Dict{String, LongDNASeq}()
        for header in headers
            ref = let
                ref = get(ref_dict, header, nothing)
                if ref === nothing
                    error("In reference JSON, header is missing: \"$(header)\"")
                else
                    ref
                end
            end
            result[ref.name] = ref.seq
            write(writer, FASTA.Record(ref.name, ref.seq))
        end
        return result
    end
end

function main(
    sample_name::AbstractString,
    cat_path::AbstractString,
    ref_dict::Dict{String, Reference},
    res_path::AbstractString,
)::Dict{String, LongDNASeq} # map of name => seq of all templates
    templates = get_templates(res_path, ref_dict)
    if isempty(templates)
        error(
            "Sample \"$sample_name\" maps to no references.
            Several steps in the pipeline assume a nonempty reference set for each sample.
            To continue, remove the following samples, delete the output, and rerun the pipeline."
        )
    end
    return dump_sequences(cat_path, templates, ref_dict)
end

end # module

"""
    IterativeAssembly

Module to iteratively assemble sequences based on the KMA aligner/assembler.

Possible todo: If you need to optimize this for speed, you can implement some
logic that takes out "completed" segments, i.e. those with 100% identity to
their template, from the loop. This will complicate matters and probably not
be worth it, though.
"""
module IterativeAssembly

using Serialization: deserialize
using FASTX
using Influenza: DEFAULT_DNA_ALN_MODEL, alignment_identity, Segment, Reference
using BioSequences: LongDNASeq, each, DNAMer, DNA, isgap
using BioAlignments: pairalign, alignment, OverlapAlignment, PairwiseAlignment
using KMATools: parse_res, parse_mat
using CodecZlib: GzipDecompressorStream

const MAX_ITERS = 6

struct Assembly
    name::String
    segment::Segment
    seq::LongDNASeq
    identity::Float64
    coverage::Float64
    depth::Float64
    depths::Vector{UInt32}
end

function main(
    sample_name::String,
    ref_dict::Dict{String, Reference},
    read_paths::Union{AbstractString, NTuple{2, AbstractString}},
    template_path::AbstractString,
    old_seq_dict::Dict{String, LongDNASeq},
    outdir::AbstractString,
    logdir::AbstractString,
    k::Int,
    convergence_threshold::AbstractFloat # should be lower for Nanopore?
)
    local mapbase
    local convergence
    iter = 0

    # Main loop - the "iterative" in "iterative assembly"
    while true
        iter += 1
        # Index template_path
        println(stderr, "Indexing iteration $iter...")
        indexbase = joinpath(outdir, "kmaindex_$(iter)")
        indexlog = joinpath(logdir, "kmaindex_$(iter)_$(sample_name).log")
        indexk = iter == 1 ? 10 : k
        index(template_path, indexbase, indexk, indexlog)

        # Map reads to index to create kma_$iter files
        println(stderr, "Mapping iteration $iter...")
        mapbase = joinpath(outdir, "kma_$(iter)")
        asm_path = mapbase * ".fsa"
        res_path = mapbase * ".res"
        mat_path = mapbase * ".mat.gz"
        maplog = joinpath(logdir, "kma_$(iter)_$(sample_name).log")
        if read_paths isa AbstractString # should be statically resolved
            kma_nanopore(read_paths, mapbase, indexbase, maplog)
        else
            kma_illumina(read_paths..., mapbase, indexbase, maplog)
        end

        # Deduplicate input
        println(stderr, "Deduplicating iteration $iter...")
        (assemblies, dedup_segments) = deduplicate(parse_fna(asm_path, res_path, mat_path, ref_dict))
        template_path = joinpath(outdir, "template_$(iter).fna")
        old_seq_dict = save_deduplicated(template_path, assemblies, old_seq_dict, dedup_segments)
        println(stderr, "Existing:", "\n\t", join(("$(a.name) $(a.identity)" for a in assemblies), "\n\t"))

        # Break if possible
        is_converged, convergence = has_converged(assemblies, dedup_segments, convergence_threshold)
        if is_converged || iter ≥ MAX_ITERS
            break
        end
    end

    # Write convergence report
    open(joinpath(outdir, "convergence.tsv"), "w") do io
        println(io, "template\tsegment\tisconverged\tid")
        for i in convergence
            println(io,
                i.template, '\t',
                i.segment, '\t',
                i.isconverged, '\t',
                round(i.identity, digits=4)
            )
        end
    end

    # Finally, rename the last one to "final"
    for ext in (".aln", ".mat.gz", ".res")
        cp(mapbase * ext, joinpath(outdir, "kma_final$(ext)"), force=true)
    end
    cp(template_path, joinpath(outdir, "kma_final.fsa"), force=true)
    cleanup(outdir, iter)
    return nothing
end

function index(
    input::AbstractString,
    outbase::AbstractString,
    k::Int,
    log::AbstractString
)
    run(pipeline(`kma index -nbp -k $k -i $input -o $outbase`, stderr=log))
end

function kma_illumina(
    infw::AbstractString,
    inrv::AbstractString,
    outbase::AbstractString,
    db::AbstractString,
    log::AbstractString
)
    cmd = `kma -ipe $infw $inrv -o $outbase -t_db $db
    -t $(Threads.nthreads()) -1t1 -mrs 0.3 -gapopen -5 -ConClave 2 -nf -matrix`
    run(pipeline(cmd, stderr=log))
    nothing
end

function kma_nanopore(
    inpath::AbstractString,
    outbase::AbstractString,
    db::AbstractString,
    log::AbstractString
)
    cmd = `kma -i $inpath -o $outbase -t_db $db
    -mp 20 -bc 0.7 -t $(Threads.nthreads()) -1t1 -mrs 0.3 -ConClave 2 -bcNano
    -nf -matrix`
    run(pipeline(cmd, stderr=log))
end

"Return vec of deduplicated, set of segments with asms removed"
function deduplicate(asms::Vector{Assembly})::Tuple{Vector{Assembly}, Set{Segment}}
    # Aligning sequences is the speed-limiting factor. We avoid aligning sequences
    # from different segments, since they should realistically never become identical.
    bysegment = Dict{Segment, Vector{Assembly}}()
    for asm in asms
        push!(get!(valtype(bysegment), bysegment, asm.segment), asm)
    end

    # If fewer than half of segments are duplicated, we have extra strict criteria
    # and are more likely to deduplicate the other segments
    few_duplicated = 2 * sum(i -> length(i) > 1, values(bysegment), init=0) < length(bysegment)

    toremove = Set{Assembly}()
    deduplicated = Dict(k => Set{Assembly}() for k in keys(bysegment))
    for asm in asms
        set = deduplicated[asm.segment]
        addcurrent = true
        empty!(toremove)
        for other in set
            bestasm = better(asm, other, few_duplicated)
            if bestasm === nothing
                nothing
            elseif bestasm === asm
                push!(toremove, other)
            elseif bestasm === other
                addcurrent = false
                break
            else
                @assert false "Unreachable"
            end
        end
        addcurrent && push!(set, asm)
        setdiff!(set, toremove)
    end
    final_segments = reduce(values(deduplicated), init=Assembly[]) do v, set
        append!(v, set)
    end
    dedup_seg = Set(k for (k,v) in deduplicated if length(v) < length(bysegment[k]))
    return (final_segments, dedup_seg)
end

# nothing if neither is better, else return the better of a or b
function better(a::Assembly, b::Assembly, few_duplicated::Bool)::Union{Nothing, Assembly}
    # These parameters are sort of arbitrary
    max_depth_ratio = few_duplicated ? 25 : 50
    max_identity = few_duplicated ? 0.97 : 0.98

    # We return the naively_better in many cases below
    naively_better = let
        if min(a.coverage, b.coverage) > 0.95
            a.depth > b.depth ? a : b
        else
            a.coverage > b.coverage ? a : b
        end
    end

    # Shortest segment is > 850 bp, so 500 is absolute minimum to ensure the assembler
    # doesn't just assemble based on a few short streches
    if min(length(a.seq), length(a.seq)) < 500
        return length(a.seq) < length(b.seq) ? b : a
    end

    # If the more common has > 50x the depth of the minor variant, I don't trust it
    # maybe this criteria is silly, but I think a 50x more abundant variant will
    # completely dominate the smaller virus in the sense that there may be as many
    # mismapped reads from the more abundant one as properly low-abundance reads
    bigdepth, smalldepth = minmax(a.depth, b.depth)
    if bigdepth / smalldepth > max_depth_ratio
        return a.depth > b.depth ? a : b
    end

    # Do a fast check for high kmer jaccard sim to avoid having to calculate
    # identity - see below for explanation
    ksim = kmer_jaccard_sim(a, b)
    if ksim > 0.75
        return naively_better
    end

    aln = alignment(pairalign(OverlapAlignment(), a.seq, b.seq, DEFAULT_DNA_ALN_MODEL))
    id = alignment_identity(OverlapAlignment(), aln)

    # This can happen if they are just too dissimilar
    if id === nothing
        return nothing
    end

    # Different part of the sequences can adapt to different references, i.e. the 3'
    # adapt to ref A and 5' to ref B. In that case the assemblies appear distinct, because
    # the parts poorly covered by reads will be different due to uncertain asssemblies in
    # those areas.
    # We check if there are BOTH big parts of the segment where A has higher depth than B, and
    # parts where B >> A, and remove a segment if so.
    if is_mutually_much_deeper(aln, a, b)
        return naively_better
    end

    # Kmer jaccard sim and sequence id follow each other with some variation. If mismatches
    # are clustered in one region of the sequences, the ksim will be disproportionately higher.
    # The thresholds below are generous thresholds based on measuring 16,000 reference pairs.
    # If we have too many assemblies at this point, the assemblies might begin to converge
    # in local part of the sequence only. This will cause ksim to shoot up and will be detected here.
    max_ksim = id > 0.95 ? 0.75 :
        id > 0.9 ? 0.5 :
        id > 0.85 ? 0.25 :
        0.2

    if ksim > max_ksim
        return naively_better
    end

    # Else if they are too different, neither one is best
    if id < max_identity
        return nothing
    end

    # Else, we simply pick the most abundant one since that will probably win
    # out the consensus sequence generation anyway through majority vote
    return naively_better
end

# Is the depths such that there are large sections where both query is much deeper
# than reference and sections where the oppotiste is true?
function is_mutually_much_deeper(aln::PairwiseAlignment, query::Assembly, ref::Assembly)::Bool
    n_q_deeper = n_r_deeper = 0
    q_index = r_index = 0
    n_matches = 0
    for (q, r) in aln
        if !isgap(q)
            q_index += 1
            # Safety measure because I'm not 100% sure when KMA assigns bases
            # as opposed to gaps.
            q_index > length(query.depths) && break
        end
        if !isgap(r)
            r_index += 1
            r_index > length(ref.depths) && break
        end
        if !isgap(r) && !isgap(q)
            n_matches += 1
            qd = query.depths[q_index]
            rd = ref.depths[r_index]
            n_q_deeper += is_much_deeper(qd, rd)
            n_r_deeper += is_much_deeper(rd, qd)
        end
    end
    min_frac_deeper, max_frac_deeper = minmax(n_q_deeper, n_r_deeper) ./ n_matches
    min_frac_deeper > 0.05 && max_frac_deeper > 0.1
end

# Whether a is significantly deeper than "than"
is_much_deeper(a::Integer, than::Integer)::Bool = (max(a, than) > 50) & (a > 5*than)

function kmer_jaccard_sim(a::Assembly, b::Assembly)
    s1 = Set(i.fw for i in each(DNAMer{31}, a.seq))
    s2 = Set(i.fw for i in each(DNAMer{31}, b.seq))
    T = min(length(s1), length(s2))
    intersect!(s1, s2)
    return length(s1) / T
end

function parse_fna(
    fsa_path::AbstractString,
    res_path::AbstractString,
    mat_path::AbstractString,
    ref_dict::Dict{String, Reference}
)::Vector{Assembly}
    res = open(io -> parse_res(io, res_path), res_path)
    res_by_header = Dict(i.template => i for i in res)
    depths_by_header = Dict(
        name => parse_depth_vector(v)
        for (name, v) in open(io -> parse_mat(io, mat_path), GzipDecompressorStream, mat_path)
    )
    open(FASTA.Reader, fsa_path) do reader
        asms = Assembly[]
        record = FASTA.Record()
        while !eof(reader)
            read!(reader, record)
            header = FASTA.header(record)::AbstractString
            row = res_by_header[header]
            depths = depths_by_header[header]
            name = String(strip(FASTA.header(record)))
            segment = ref_dict[name].segment
            seq = FASTA.sequence(LongDNASeq, record)
            
            # Query and template identity are calculated differently, because gap
            # positions in query/template are not counted towards identity.
            # Hence an upper bound of the total alignment identity is the minimum
            # of the two ids
            id = min(row.tid, row.qid)
            push!(asms, Assembly(name, segment, seq, id, row.tcov, row.depth, depths))
        end
        return asms
    end
end

function parse_depth_vector(depths_vector::Vector{<:Tuple{DNA, Tuple{Vararg{Integer}}}})
    res = UInt32[]
    for (_dna, deps) in depths_vector
        n_gaps = last(deps)
        biggest = maximum(deps[1:lastindex(deps)-1])
        if n_gaps ≤ biggest
            nongaps = sum(deps[1:lastindex(deps)-1], init=zero(eltype(deps)))
            push!(res, UInt32(nongaps))
        end
    end
    res
end

function has_converged(
    assemblies::Vector{Assembly},
    dedup_segment::Set{Segment},
    convergence_threshold::AbstractFloat
)::Tuple{Bool, Vector{<:NamedTuple{(:template, :segment, :isconverged, :identity)}}}
    convergence = map(assemblies) do assembly
        template = assembly.name
        isconverged = let
            assembly.identity ≥ convergence_threshold &&
            # Non-1.0 coverage implies indels, which often linger around and necessitate
            # additional rounds. We allow float rounding error here.
            assembly.coverage > 0.999999 &&
            assembly.coverage < 1.000001
        end
        identity = assembly.identity
        segment = assembly.segment
        (; template, segment, isconverged, identity)
    end
    isconverged = isempty(dedup_segment) && all(i -> i.isconverged, convergence)
    return (isconverged, convergence)
end

# If a segment was deduplicated, we keep using the old assembly instead of replacing
# it with the new assembly. This is because, when there are redundant sequences present,
# the assemblies are often of poor quality as the reads are spread too thin.
function save_deduplicated(
    outpath::AbstractString,
    assemblies::Vector{Assembly},
    old_seq_dict::Dict{String, LongDNASeq},
    dedup_segments::Set{Segment}
)::Dict{String, LongDNASeq}
    open(FASTA.Writer, outpath) do writer
        result = Dict{String, LongDNASeq}()
        for assembly in assemblies
            seq = if assembly.segment in dedup_segments
                old_seq_dict[assembly.name]
            else
                assembly.seq
            end
            write(writer, FASTA.Record(assembly.name, seq))
            result[assembly.name] = seq
        end
        return result
    end
end

function cleanup(dir::AbstractString, iters::Integer)
    # We don't delete the first iteration.
    # The final iteration has been copied to `final`
    for i in 1:iters
        i > 1 && rm(joinpath(dir, "kma_$(i).aln"))
        i > 1 && rm(joinpath(dir, "kma_$(i).mat.gz"))
        i > 1 && rm(joinpath(dir, "kma_$(i).fsa"))
        i > 1 && rm(joinpath(dir, "template_$(i-1).fna"))
        rm(joinpath(dir, "kmaindex_$(i).comp.b"))
        rm(joinpath(dir, "kmaindex_$(i).length.b"))
        rm(joinpath(dir, "kmaindex_$(i).name"))
        rm(joinpath(dir, "kmaindex_$(i).seq.b"))
    end
    rm(joinpath(dir, "template_initial.fna"))
end

end # module

using Influenza: Influenza

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) ∉ (9, 10)
        error("Usage: julia iter_asm.jl samplename refjson templatepath respath outdir logdir k threshold read1 [read2]")
    end
    samplename, jsonpath, templatepath, respath, outdir, logdir = ARGS[1:6]
    refdict = let
        v = Influenza.load_references(jsonpath)
        Dict(i.name => i for i in v)
    end
    k = parse(Int, ARGS[7])
    threshold = parse(Float64, ARGS[8])
    if threshold < 0.0 || threshold > 1.0
        error("Threshold must be in 0.0-1.0")
    end
    readpaths = length(ARGS) == 9 ? ARGS[9] : (ARGS[9], ARGS[10])

    old_seq_dict = SelectTemplates.main(samplename, templatepath, refdict, respath)
    IterativeAssembly.main(
        samplename, refdict, readpaths, templatepath, old_seq_dict, outdir, logdir, k, threshold
    )
end