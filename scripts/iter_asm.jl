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
using Influenza: DEFAULT_DNA_ALN_MODEL, alignment_identity, Segment
using FluWorkflowTools: split_segment
using BioSequences: LongDNASeq
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

FASTA.Record(a::Assembly) = FASTA.Record("$(a.name)_$(a.segment)", a.seq)

function main(
    samplename::String,
    readpaths::Union{AbstractString, NTuple{2, AbstractString}},
    asm_path::AbstractString,
    res_path::AbstractString,
    outdir::AbstractString,
    logdir::AbstractString,
    k::Int,
    convergence_threshold::AbstractFloat # should be lower for Nanopore?
)
    # First iteration has been run outside this script with slightly different params
    iter = 1
    dedup_path = joinpath(outdir, "deduplicated_$(iter).fna")
    assemblies = deduplicate_and_save(parse_fna(asm_path, res_path), dedup_path)
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
        assemblies = deduplicate_and_save(parse_fna(asm_path, res_path), dedup_path)
        println(stderr, "Existing:", "\n\t", join(("$(a.name) $(a.identity)" for a in assemblies), "\n\t"))

        # Break if possible
        if has_conveged(assemblies, convergence_threshold) || iter == MAX_ITERS
            break
        end

        iter += 1
    end

    # Finally, rename the last one to "final"
    for ext in (".aln", ".fsa", ".mat.gz", ".res")
        mv(joinpath(outdir, "kma_$(iter)$(ext)"), joinpath(outdir, "kma_final$(ext)"))
    end

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
    -t $(Threads.nthreads()) -1t1 -gapopen -5 -nf -matrix`
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
    -mp 20 -bc 0.7 -t $(Threads.nthreads()) -1t1 -bcNano -nf -matrix`
    run(pipeline(cmd, stderr=log))
end

function deduplicate_and_save(asms::Vector{Assembly}, path::AbstractString)::Vector{Assembly}
    next = deduplicate(asms)
    open(FASTA.Writer, path) do writer
        foreach(a -> write(writer, FASTA.Record(a)), asms)
    end
    return next
end

function deduplicate(asms::Vector{Assembly})::Vector{Assembly}
    # Aligning sequences is the speed-limiting factor. We avoid aligning sequences
    # from different segments, since they should realistically never become identical.
    bysegment = Dict{Segment, Vector{Assembly}}()
    for asm in asms
        push!(get!(valtype(bysegment), bysegment, asm.segment), asm)
    end

    toremove = Set{Assembly}()
    deduplicated = Dict(k => Set{Assembly}() for k in keys(bysegment))
    for asm in asms
        set = deduplicated[asm.segment]
        addcurrent = true
        empty!(toremove)
        for other in set
            bestasm = better(asm, other)
            if bestasm === nothing
                nothing
            elseif bestasm === asm
                push!(toremove, other)
            elseif bestasm === other
                addcurrent = false
                break
            else
                error() # unreachable!
            end
        end
        addcurrent && push!(set, asm)
        setdiff!(set, toremove)
    end
    return reduce(values(deduplicated), init=Assembly[]) do v, set
        append!(v, set)
    end
end

function better(a::Assembly, b::Assembly)::Union{Nothing, Assembly}
    # Shortest segment is > 850 bp, so 800 is absolute minimum to ensure the assembler
    # doesn't just assemble based on a few short streches
    if min(length(a.seq), length(a.seq)) < 800
        return length(a.seq) < length(b.seq) ? b : a
    end

    aln = alignment(pairalign(OverlapAlignment(), a.seq, b.seq, DEFAULT_DNA_ALN_MODEL))
    id = alignment_identity(OverlapAlignment(), aln)
    if isnothing(id) || id < 0.98 # this is sort of arbitrary, maybe it should be a parameter
        return nothing
    end

    return a.depth < b.depth ? b : a
end

function parse_fna(fsa_path::AbstractString, res_path::AbstractString)::Vector{Assembly}
    res = open(io -> parse_res(io, res_path), res_path)
    res_by_header = Dict(i.template => i for i in res)
    open(FASTA.Reader, fsa_path) do reader
        asms = Assembly[]
        record = FASTA.Record()
        while !eof(reader)
            read!(reader, record)
            header = FASTA.header(record)
            row = res_by_header[header]
            name, segment = split_segment(FASTA.header(record))
            seq = FASTA.sequence(record)
            push!(asms, Assembly(name, segment, seq, row.tid, row.tcov, row.depth))
        end
        return asms
    end
end

function has_conveged(assemblies::Vector{Assembly}, convergence_threshold)::Bool
    all(assemblies) do assembly
        assembly.identity ≥ convergence_threshold
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) ∉ (8, 9)
        error("Usage: julia iter_asm.jl samplename asmpath respath outdir logdir k threshold read1 [read2]")
    end
    samplename, asmpath, respath, outdir, logdir = ARGS[1:5]
    k = parse(Int, ARGS[6])
    threshold = parse(Float64, ARGS[7])
    if threshold < 0.0 || threshold > 1.0
        error("Threshold must be in 0.0-1.0")
    end
    readpaths = length(ARGS) == 8 ? ARGS[8] : (ARGS[8], ARGS[9])
    main(
        samplename,
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

