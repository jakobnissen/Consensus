"""
    iter_asm

Script to iteratively assemble sequences based on the KMA aligner/assembler.

Possible todo: If you need to optimize this for speed, you can implement some
logic that takes out "completed" segments, i.e. those with 100% identity to
their template, from the loop. This will complicate matters and probably not
be worth it, though.
"""
module t

using FASTX
using Influenza: DEFAULT_DNA_ALN_MODEL, alignment_identity, Segment, load_references
using BioSequences: LongDNASeq, each, DNAMer
using BioAlignments: pairalign, alignment, OverlapAlignment
using KMATools: parse_res

const MAX_ITERS = 6

struct Assembly
    name::String
    segment::Segment
    seq::LongDNASeq
    identity::Float64
    coverage::Float64
    depth::Float64
end

function main(
    samplename::String,
    jsonpath::AbstractString,
    readpaths::Union{AbstractString, NTuple{2, AbstractString}},
    asm_path::AbstractString,
    res_path::AbstractString,
    outdir::AbstractString,
    logdir::AbstractString,
    k::Int,
    convergence_threshold::AbstractFloat # should be lower for Nanopore?
)
    segment_map = Dict(ref.name => ref.segment for ref in load_references(jsonpath))
    # First iteration has been run outside this script with slightly different params
    iter = 1
    dedup_path = joinpath(outdir, "deduplicated_$(iter).fna")
    (assemblies, has_deduplicated) = deduplicate_and_save(parse_fna(asm_path, res_path, segment_map), dedup_path)
    local mapbase
    local convergence
    iter += 1
    while true
        # Index
        println(stderr, "Indexing iteration $iter...")
        indexbase = joinpath(outdir, "kmaindex_$(iter)")
        indexlog = joinpath(logdir, "kmaindex_$(iter)_$(samplename).log")
        index(dedup_path, indexbase, k, indexlog)

        # Map
        println(stderr, "Mapping iteration $iter...")
        mapbase = joinpath(outdir, "kma_$(iter)")
        asm_path = mapbase * ".fsa"
        res_path = mapbase * ".res"
        maplog = joinpath(logdir, "kma_$(iter)_$(samplename).log")
        if readpaths isa String # should be statically resolved
            kma_nanopore(readpaths, mapbase, indexbase, maplog)
        else
            kma_illumina(readpaths..., mapbase, indexbase, maplog)
        end

        # Deduplicate input
        println(stderr, "Deduplicating iteration $iter...")
        dedup_path = joinpath(outdir, "deduplicated_$(iter).fna")
        (assemblies, has_deduplicated) = deduplicate_and_save(parse_fna(asm_path, res_path, segment_map), dedup_path)
        println(stderr, "Existing:", "\n\t", join(("$(a.name) $(a.identity)" for a in assemblies), "\n\t"))

        # Break if possible
        is_converged, convergence = has_conveged(assemblies, has_deduplicated, convergence_threshold)
        if is_converged || iter ≥ MAX_ITERS
            break
        end

        iter += 1
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
    cp(dedup_path, joinpath(outdir, "kma_final.fsa"), force=true)
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
    -t $(Threads.nthreads()) -1t1 -mrs 0.3 -gapopen -5 -nf -matrix`
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
    -mp 20 -bc 0.7 -t $(Threads.nthreads()) -1t1 -mrs 0.3 -bcNano -nf -matrix`
    run(pipeline(cmd, stderr=log))
end

function deduplicate_and_save(
    asms::Vector{Assembly},
    path::AbstractString
)::Tuple{Vector{Assembly}, Bool}
    next = deduplicate(asms)
    open(FASTA.Writer, path) do writer
        foreach(a -> write(writer, FASTA.Record(a.name, a.seq)), next)
    end
    has_deduplicated = length(next) != length(asms)
    return (next, has_deduplicated)
end

function deduplicate(asms::Vector{Assembly})::Vector{Assembly}
    # Aligning sequences is the speed-limiting factor. We avoid aligning sequences
    # from different segments, since they should realistically never become identical.
    bysegment = Dict{Segment, Vector{Assembly}}()
    for asm in asms
        push!(get!(valtype(bysegment), bysegment, asm.segment), asm)
    end

    # If fewer than half of segments are duplicated, we have extra strict criteria
    # and are more likely to deduplicate the other segments
    few_duplicated = 2 * sum(i -> length(i) > 1, values(bysegment)) < length(bysegment)

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
    return reduce(values(deduplicated), init=Assembly[]) do v, set
        append!(v, set)
    end
end

function better(a::Assembly, b::Assembly, few_duplicated::Bool)::Union{Nothing, Assembly}
    # These parameters are sort of arbitrary
    max_depth_ratio = few_duplicated ? 25 : 100
    max_identity = few_duplicated ? 0.97 : 0.98

    naively_better = let
        if min(a.coverage, b.coverage) > 0.98
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

    # If the more common has > 100x the depth of the minor variant, I don't trust it
    # maybe this criteria is silly, but I think a 100x more abundant variant will
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

    # This can happen if one of the sequences are absolutely horribly reconstructed
    if id === nothing
        return naively_better
    end

    # Kmer jaccard sim and sequence id follow each other with some variation. If mismatches
    # are clustered in one region of the sequences, the ksim will be disproportionately higher.
    # The thresholds below are generous thresholds based on checking 16,000 references.
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
    if id === nothing || id < max_identity
        return nothing
    end

    # Else, we simply pick the most abundant one since that will probably win
    # out the consensus sequence generation anyway through majority vote
    return naively_better
end

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
    segment_map::Dict{String, Segment}
)::Vector{Assembly}
    res = open(io -> parse_res(io, res_path), res_path)
    res_by_header = Dict(i.template => i for i in res)
    open(FASTA.Reader, fsa_path) do reader
        asms = Assembly[]
        record = FASTA.Record()
        while !eof(reader)
            read!(reader, record)
            header = FASTA.header(record)
            row = res_by_header[header]
            name = String(strip(FASTA.header(record)))
            segment = segment_map[name]
            seq = FASTA.sequence(LongDNASeq, record)
            
            # Query and template identity are calculated differently, because gap
            # positions in query/template are not counted towards identity.
            # Hence an upper bound of the total alignment identity is the minimum
            # of the two ids
            id = min(row.tid, row.qid)
            push!(asms, Assembly(name, segment, seq, id, row.tcov, row.depth))
        end
        return asms
    end
end

function has_conveged(
    assemblies::Vector{Assembly},
    has_deduplicated::Bool,
    convergence_threshold::AbstractFloat
)
    convergence = map(assemblies) do assembly
        template = assembly.name
        isconverged = assembly.identity ≥ convergence_threshold
        identity = assembly.identity
        segment = assembly.segment
        (; template, segment, isconverged, identity)
    end
    isconverged = !has_deduplicated && all(i -> i.isconverged, convergence)
    return (isconverged, convergence)
end

function cleanup(dir::AbstractString, iters::Integer)
    for i in 2:iters
        rm(joinpath(dir, "kma_$(i).aln"))
        rm(joinpath(dir, "kma_$(i).fsa"))
        rm(joinpath(dir, "kma_$(i).mat.gz"))
        rm(joinpath(dir, "kmaindex_$(i).comp.b"))
        rm(joinpath(dir, "kmaindex_$(i).length.b"))
        rm(joinpath(dir, "kmaindex_$(i).name"))
        rm(joinpath(dir, "kmaindex_$(i).seq.b"))
        rm(joinpath(dir, "deduplicated_$(i-1).fna"))
    end
    rm(joinpath(dir, "kma_$(iters).res"))
end


if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) ∉ (9, 10)
        error("Usage: julia iter_asm.jl samplename jsonpath asmpath respath outdir logdir k threshold read1 [read2]")
    end
    samplename, jsonpath, asmpath, respath, outdir, logdir = ARGS[1:6]
    k = parse(Int, ARGS[7])
    threshold = parse(Float64, ARGS[8])
    if threshold < 0.0 || threshold > 1.0
        error("Threshold must be in 0.0-1.0")
    end
    readpaths = length(ARGS) == 9 ? ARGS[9] : (ARGS[10], ARGS[10])
    main(
        samplename,
        jsonpath,
        readpaths,
        asmpath,
        respath,
        outdir,
        logdir,
        k,
        threshold
    )
end

end # module