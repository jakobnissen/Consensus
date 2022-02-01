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
        for i in eachindex(samples)
            v = report(io, samples[i], order[i], alnasms[i], depths[i], read_stats[i], is_illumina)
            push!(result, v)
        end
        print(io, '\n')
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
