const KMerSet = Set{DNAMer{16}}

function snakemake_entrypoint(
    report_path::AbstractString, # output report
    ref_dir::AbstractString, # dir of .fna + .json ref files
    aln_dir::AbstractString, # dir of kma aln
    seq_dir::AbstractString,
    tmp_dir::AbstractString,
    is_illumina::Bool,
)::Nothing
    samples = map(Sample, sort!(readdir(aln_dir)))
    paths = map(samples) do sample
        name = nameof(sample)
        (;
            sample=sample,
            fastp=joinpath(tmp_dir, "trim", name, "report.json"),
            asm=joinpath(aln_dir, name, "kma_final.fsa"),
            convergence=joinpath(aln_dir, name, "convergence.tsv"),
            t_depth=joinpath(aln_dir, name, "kma_1.mat.gz"),
            a_depth=joinpath(aln_dir, name, "kma_final.mat.gz"),
            seq_dir=joinpath(seq_dir, name)
        )
    end
    read_stats = map(p -> check_reads(p.fastp), paths)
    aln_asms = load_aligned_assemblies(
        [s.asm for s in paths],
        joinpath(ref_dir, "refs.json")
    )
    foreach(zip(aln_asms, paths)) do (alnasmv, path)
        convergence_check(alnasmv, path.convergence)
    end
    depths = map(zip(aln_asms, paths)) do (alnasmv, path)
        load_depths_and_errors(alnasmv, path.t_depth, path.a_depth)
    end
    segmentorder = map(i -> order_alnasms(i...), zip(aln_asms, depths))

    # Re-order the vectors based on segment and order.
    for i in eachindex(aln_asms, depths, segmentorder)
        ord = sortperm([(a.reference.segment, o) for (a, o) in zip(aln_asms[i], segmentorder[i])])
        aln_asms[i] = aln_asms[i][ord]
        depths[i] = depths[i][ord]
        segmentorder[i] = segmentorder[i][ord]
    end

    passes = report(report_path, samples, segmentorder, aln_asms, depths, read_stats, is_illumina)

    for (p, alnasmv, pass, sorder) in zip(paths, aln_asms, passes, segmentorder)
        write_sequences(p.seq_dir, p.sample, alnasmv, pass, sorder)
    end

    open(GzipCompressorStream, joinpath(tmp_dir, "internal.jls.gz"), "w") do io
        serialize(io, samples, aln_asms, depths, passes, segmentorder)
    end

    return nothing
end

# This function has a bit weird logical flow, because it's easiest to print
# at the same time as determining whether the segment has passed.
# But we also need to prefix failed segments with "FAIL".
# So what we do is to print to an IOBuffer, then prefix any segment with FAIL.
function report(
    report_path::AbstractString,
    samples::Vector{Sample},
    order::Vector{Vector{UInt8}},
    alnasms::Vector{Vector{AlignedAssembly}},
    depths::Vector{Vector{Depths}},
    read_stats::Vector{<:ReadStats},
    is_illumina::Bool,
)::Vector{Vector{Bool}}
    result = Vector{Bool}[]
    open(report_path, "w") do io
        # Print report for each segment
        for i in eachindex(samples)
            v = report(io, samples[i], order[i], alnasms[i], depths[i], read_stats[i], is_illumina)
            push!(result, v)
        end
        print(io, '\n')

        # Check for duplicate segments
        check_duplicates(io, samples, order, alnasms)
    end
    return result
end

function report(
    io::IO,
    sample::Sample,
    order::Vector{UInt8},
    alnasms::Vector{AlignedAssembly},
    depths::Vector{Depths},
    read_stats::ReadStats,
    is_illumina::Bool,
)::Vector{Bool}
    indexof = Dict(a => i for (i, a) in enumerate(alnasms))
    passes = fill(false, length(alnasms))

    primary = Vector{Union{Tuple{AlignedAssembly, Depths}, Nothing}}(nothing, N_SEGMENTS)
    aux = Tuple{AlignedAssembly, Depths, UInt8}[]
    for i in eachindex(alnasms, depths, order)
        a, d, o = alnasms[i], depths[i], order[i]
        segment = a.reference.segment
        index = Integer(segment) + 0x01
        if o == 1
            primary[index] = (a, d)
        else
            push!(aux, (a, d, o))
        end
    end
    sort!(aux, by=i -> first(i).reference.segment)

    # Sample header
    println(io, sample)
    report(io, read_stats)

    # Primary segments report
    for (i, data) in enumerate(primary)
        segment = Segment(i - 1)
        (buf, pass) = report_segment(data, is_illumina)
        if data !== nothing
            alnasm, _ = data
            passes[indexof[alnasm]] = pass
        end
        print(io,
            '\t', pass ? "     " : "FAIL ",
            rpad(string(segment) * ":", 4),
        )
        write(io, take!(buf))
    end

    # Secondary segments
    if !isempty(aux)
        println(io, "")
        println(io, "\t[ POSSIBLE SUPERINFECTION ]")
        for i in aux
            alnasm, depth, order = i
            _, pass = report_segment((alnasm, depth), is_illumina)
            segment = alnasm.reference.segment
            print(io,
                "\t", pass ? "     " : "FAIL ",
                rpad(string(segment), 3) * ' ' * string(order) * ':'
            )
            print_segment_header(io, alnasm, depth)
            passes[indexof[alnasm]] = pass
        end
    end

    print(io, '\n')
    return passes
end

function report_segment(
    data::Union{Tuple{AlignedAssembly, Depths}, Nothing}, is_illumina::Bool
)::Tuple{IOBuffer, Bool}
    buf = IOBuffer()
    passed = true
    if data === nothing
        println(buf, " Missing segment")
        return (buf, false)
    end
    (alnasm, depth) = data
    print_segment_header(buf, alnasm, depth)

    # If an error is terminating, it's so bad that we don't need to print
    # the other errors as they would just be spam. E.g. if the depth is < 5
    # or coverage is < 0.75
    critical_errors = filter(alnasm.errors) do error
        (error isa Influenza.ErrorLowDepthBases && error.n > 500) ||
        (error isa Influenza.ErrorLowCoverage && error.coverage < 0.75)
    end
    if !isempty(critical_errors)
        println(buf, "\t\tERROR ", first(critical_errors))
        println(buf, "\t\t      Skipping additional errors...")
        return (buf, false)
    end

    for error in alnasm.errors
        p = pass(error, is_illumina)
        println(buf, "\t\t", (p ? "      " : "ERROR "), error)
        passed &= p
    end

    IndelError = Union{Influenza.ErrorFrameShift, Influenza.ErrorIndelTooBig}

    for protein in alnasm.proteins
        # There can be a ton of indel errors. We're not interested in printing all of them,
        # so if there are a lot of indel errors, we just replace it with a ErrorTooManyIndels.
        if count(i -> i isa IndelError, protein.errors) > 3
            N = length(protein.errors)
            filter!(i -> !isa(i, IndelError), protein.errors)
            n_indel_errors = N - length(protein.errors)
            push!(protein.errors, ErrorTooManyIndels(n_indel_errors))
        end

        for error in protein.errors
            p = pass(error, is_illumina)
            println(buf, "\t\t", (p ? "      " : "ERROR "), protein.variant, ": ", error)
            passed &= (p || !is_important(protein.variant))
        end
    end
    return (buf, passed)
end

function print_segment_header(io::IO, alnasm::AlignedAssembly, depth::Depths)
    print(io, " Identity $(@sprintf "%.3f" (alnasm.identity)),")
    println(io, " depth $(@sprintf "%.2e" assembly_depth(depth)), coverage $(@sprintf "%.3f" template_coverage(depth))")
end

# fallback: segment errors fails
pass(::Influenza.InfluenzaError, is_illumina::Bool) = false

# Seems weird, but NA stalk can have 75 bp without any problems
pass(x::Influenza.ErrorIndelTooBig, _::Bool) = length(x.indel) < 100
pass(x::Influenza.ErrorEarlyStop, _::Bool) = x.observed_naa + 14 > x.expected_naa
pass(x::Influenza.ErrorInsignificant, is_illumina::Bool) = x.n_insignificant < ifelse(is_illumina, 5, 25)
pass(x::Influenza.ErrorAmbiguous, is_illumina::Bool) = x.n_ambiguous < ifelse(is_illumina, 5, 25)
pass(x::Influenza.ErrorLateStop, is_illumina::Bool) = true

# We fail with low depth because cross-contamination from other samples
# can give enough reads to reconstruct a segment.
# If more than 500 bp have lower than 25 depth, we can't trust it.
pass(x::Influenza.ErrorLowDepthBases, is_illumina::Bool) = x.n < 500

function check_duplicates(
    io::IO,
    samples::Vector{Sample},
    order::Vector{Vector{UInt8}},
    alnasms::Vector{Vector{AlignedAssembly}},
)
    pairs = check_duplicates(samples, order, alnasms)
    isempty(pairs) && return nothing
    println(io, "Duplicate segments (> 99.5% identity)")
    for (segment, v) in sort!(collect(pairs), by=first)
        println(io, segment)
        for ((s1, ord1), (s2, ord2), id) in v
            println(io, "\t$s1 ($ord1)\t$s2 ($ord2)\t$(@sprintf "%.3f" id)")
        end
        println(io)
    end
end

function check_duplicates(
    samples::Vector{Sample},
    order::Vector{Vector{UInt8}},
    alnasms::Vector{Vector{AlignedAssembly}}
   # vector of: (segment, [((sample1, ord), (sample2, ord), id) ... ]
)::Dict{Segment, Vector{Tuple{Tuple{Sample, UInt8}, Tuple{Sample, UInt8}, Float64}}}
    bysegment = Dict{Segment, Vector{Tuple{Sample, UInt8, LongDNASeq, KMerSet}}}()
    for (sample, orders, va) in zip(samples, order, alnasms), (ord, alnasm) in zip(orders, va)
        seq = alnasm.assembly.seq
        tup = (sample, ord, seq, kmerset(seq))
        push!(get!(valtype(bysegment), bysegment, alnasm.reference.segment), tup)
    end

    result = Dict{Segment, Vector{Tuple{Tuple{Sample, UInt8}, Tuple{Sample, UInt8}, Float64}}}()
    for (segment, v) in bysegment
        for i in 1:length(v)-1, j in i+1:length(v)
            seqa, seqb = v[i][3], v[j][3]
            len1, len2 = minmax(length(seqa), length(seqb))
            len1 / len2 < 0.8 && continue
            overlap(v[i][4], v[j][4]) > 0.8 || continue
            id = alnid(seqa, seqb)
            id > 0.995 || continue
            element = ((v[i][1], v[i][2]), (v[j][1], v[j][2]), id)
            push!(get!(valtype(result), result, segment), element)
        end
    end
    return result
end

function kmerset(seq::LongDNASeq)::KMerSet
    Set(canonical(mer) for mer in each(DNAMer{16}, seq))
end

function overlap(a::T, b::T) where {T <: Set{<:DNAMer}}
    length(intersect(a, b)) / min(length(a), length(b))
end

function alnid(a::LongDNASeq, b::LongDNASeq)::Float64
    # This is never nothing because we assume the alignment doesnt fail
    aln = BioAlignments.alignment(
        pairalign(OverlapAlignment(), a, b, DEFAULT_DNA_ALN_MODEL)
    )::BioAlignments.PairwiseAlignment
    # And this is only nothing if the identity is very low, but we already
    # checked with kmer overlap, so this should never happen
    Influenza.alignment_identity(OverlapAlignment(), aln)::Float64
end
