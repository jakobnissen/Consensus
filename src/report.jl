function illumina_snakemake_entrypoint(
    report_path::AbstractString, # output report
    ref_dir::AbstractString, # dir of .fna + .json ref files
    aln_dir::AbstractString, # dir of medaka / kma aln
    seq_dir::AbstractString,
    depths_plot_dir::AbstractString,
    tmp_dir::AbstractString
)::Nothing
    samples = map(Sample, sort!(readdir(aln_dir)))
    paths = map(samples) do sample
        name = nameof(sample)
        (;
            sample=sample,
            fastp=joinpath(tmp_dir, "trim", name, "report.json"),
            asm=joinpath(aln_dir, name, "kma2.fsa"),
            res=joinpath(aln_dir, name, "kma2.res"),
            t_depth=joinpath(aln_dir, name, "kma1.mat.gz"),
            a_depth=joinpath(aln_dir, name, "kma2.mat.gz"),
            t_plot=joinpath(depths_plot_dir, name * "_template.pdf"),
            a_plot=joinpath(depths_plot_dir, name * "_assembly.pdf"),
            seq_dir=joinpath(seq_dir, name)
        )
    end
    read_stats = map(p -> check_reads(p.fastp), paths)
    aln_asms = load_aligned_assemblies(
        [s.asm for s in paths],
        [s.sample for s in paths],
        joinpath(ref_dir, "refs.json"),
        true
    )
    foreach(zip(aln_asms, paths)) do (alnasmv, path)
        kma2_identity_check(alnasmv, path.res)
    end
    depths = map(zip(aln_asms, paths)) do (alnasmv, path)
        load_depths_and_errors(alnasmv, path.t_depth, path.a_depth)
    end
    isprimary = map(i -> get_primary(i...), zip(aln_asms, depths))
    
    for (pr, path, de, alnasmv) in zip(isprimary, paths, depths, aln_asms)
        v = [(alnasmv[i].reference.segment, de[i]) for i in eachindex(de, alnasmv, pr) if pr[i]]
        plot_depths(path.t_plot, path.a_plot, v)
    end

    passes = report(report_path, samples, isprimary, aln_asms, depths, read_stats)

    for (p, alnasmv, pass, prim) in zip(paths, aln_asms, passes, isprimary)
        write_sequences(p.seq_dir, p.sample, alnasmv, pass, prim)
    end

    open(joinpath(tmp_dir, "internal.json"), "w") do io
        serialize_alnasms(io, aln_asms, passes, isprimary)
    end

    return nothing
end

function nanopore_snakemake_entrypoint(
    report_path::AbstractString, # output report
    ref_dir::AbstractString, # dir of .fna + .json ref files
    aln_dir::AbstractString, # dir of medaka / kma aln
    seq_dir::AbstractString,
    depths_plot_dir::AbstractString,
    tmp_dir::AbstractString
)::Nothing
    samples = map(Sample, sort!(readdir(aln_dir)))
    paths = map(samples) do sample
        name = nameof(sample)
        (;
            sample=sample,
            fastp=joinpath(tmp_dir, "trim", name, "report.json"),
            asm=joinpath(aln_dir, name, "medaka", "consensus.fasta"),
            t_depth=joinpath(aln_dir, name, "kma1.mat.gz"),
            t_plot=joinpath(depths_plot_dir, name * "_template.pdf"),
            seq_dir=joinpath(seq_dir, name)
        )
    end
    read_stats = map(p -> check_reads(p.fastp), paths)
    aln_asms = load_aligned_assemblies(
        [s.asm for s in paths],
        [s.sample for s in paths],
        joinpath(ref_dir, "refs.json"),
        false
    )

    depths = map(zip(aln_asms, paths)) do (alnasmv, path)
        load_depths_and_errors(alnasmv, path.t_depth)
    end
    isprimary = map(i -> get_primary(i...), zip(aln_asms, depths))
    
    for (pr, path, de, alnasmv) in zip(isprimary, paths, depths, aln_asms)
        v = [(alnasmv[i].reference.segment, de[i]) for i in eachindex(de, alnasmv, pr) if pr[i]]
        plot_depths(path.t_plot, path.t_plot, v)
    end

    passes = report(report_path, samples, isprimary, aln_asms, depths, read_stats)

    for (p, alnasmv, pass, prim) in zip(paths, aln_asms, passes, isprimary)
        write_sequences(p.seq_dir, p.sample, alnasmv, pass, prim)
    end

    open(joinpath(tmp_dir, "internal.json"), "w") do io
        serialize_alnasms(io, aln_asms, passes, isprimary)
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
    isprimary::Vector{Vector{Bool}},
    alnasms::Vector{Vector{AlignedAssembly}},
    depths::Vector{Vector{Depths}},
    read_stats::Vector{<:ReadStats},
)::Vector{Vector{Bool}}
    result = Vector{Bool}[]
    open(report_path, "w") do io
        for i in eachindex(samples)
            v = report(io, samples[i], isprimary[i], alnasms[i], depths[i], read_stats[i])
            push!(result, v)
        end
        print(io, '\n')
    end
    return result
end

function report(
    io::IO,
    sample::Sample,
    isprimary::Vector{Bool},
    alnasms::Vector{AlignedAssembly},
    depths::Vector{Depths},
    read_stats::ReadStats,
)::Vector{Bool}
    indexof = Dict(a => i for (i, a) in enumerate(alnasms))
    passes = fill(false, length(alnasms))

    primary = Vector{Union{Tuple{AlignedAssembly, Depths}, Nothing}}(undef, N_SEGMENTS)
    fill!(primary, nothing)
    aux = Tuple{AlignedAssembly, Depths}[]
    for i in eachindex(alnasms, depths, isprimary)
        a, d, p = alnasms[i], depths[i], isprimary[i]
        segment = a.reference.segment
        index = Integer(segment) + 0x01
        if p
            primary[index] = (a, d)
        else
            push!(aux, (a, d))
        end
    end
    sort!(aux, by=i -> first(i).reference.segment)

    # Sample header
    println(io, sample)
    report(io, read_stats)

    # Primary segments report
    for (i, data) in enumerate(primary)
        segment = Segment(i - 1)
        (buf, pass) = report_segment(data)
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
            alnasm, depth = i
            _, pass = report_segment(i)
            segment = alnasm.reference.segment
            print(io,
                "\t", pass ? "     " : "FAIL ",
                rpad(string(segment) * ":", 4)
            )
            print_segment_header(io, alnasm, depth)
            passes[indexof[alnasm]] = pass
        end
    end

    print(io, '\n')
    return passes
end

function report_segment(
    data::Union{Tuple{AlignedAssembly, Depths}, Nothing}
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
        p = pass(error)
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
            p = pass(error)
            println(buf, "\t\t", (p ? "      " : "ERROR "), protein.variant, ": ", error)
            passed &= (p || !is_important(protein.variant))
        end
    end
    return (buf, passed)
end

function print_segment_header(io::IO, alnasm::AlignedAssembly, depth::Depths)
    print(io, " Identity $(@sprintf "%.3f" (alnasm.identity)),")
    println(io, " depth $(@sprintf "%.2e" mean_depth(depth)), coverage $(@sprintf "%.3f" coverage(depth))")
end

# fallback: segment errors fails
pass(::Influenza.InfluenzaError) = false
pass(x::Influenza.ErrorEarlyStop) = x.observed_naa + 14 > x.expected_naa
pass(x::Influenza.ErrorInsignificant) = x.n_insignificant < 5
pass(x::Influenza.ErrorAmbiguous) = x.n_ambiguous < 5
pass(x::Influenza.ErrorLateStop) = true

# Currently, we can't correctly estimate depths of Nanopore reads (see depths.jl file).
# hence, to avoid getting too many false failures, we allow this error to pass through.
# Too low depth will almost certainly be reflected in other errors anyway
pass(x::Influenza.ErrorLowDepthBases) = true
