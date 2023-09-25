const KMerSet = Set{DNAMer{16}}

function snakemake_entrypoint(
    config_path::AbstractString,
    report_path::AbstractString, # output report
    ref_dir::AbstractString, # dir of .fna + .json ref files
    aln_dir::AbstractString, # dir of kma aln
    seq_dir::AbstractString,
    depths_dir::AbstractString,
    tmp_dir::AbstractString,
    similar::Bool,
    is_illumina::Bool,
)::Nothing
    config = Config(Dict(JSON3.read(config_path)), is_illumina)

    samples = map(Sample, sort!(filter!(i -> !startswith(i, '.'), readdir(aln_dir))))
    paths = map(samples) do sample
        name = nameof(sample)
        (;
            sample=sample,
            fastp=joinpath(tmp_dir, "trim", name, "report.json"),
            asm=joinpath(aln_dir, name, "kma_final.fsa"),
            convergence=joinpath(aln_dir, name, "convergence.tsv"),
            t_depth=joinpath(aln_dir, name, "kma_1.mat.gz"),
            a_depth=joinpath(aln_dir, name, "kma_final.mat.gz"),
            seq_dir=joinpath(seq_dir, name),
        )
    end
    read_stats = map(p -> check_reads(p.fastp), paths)
    aln_asms =
        load_aligned_assemblies([s.asm for s in paths], joinpath(ref_dir, "refs.json"))
    foreach(zip(aln_asms, paths)) do (alnasmv, path)
        convergence_check(alnasmv, path.convergence)
    end
    depths = Vector{Vector{Depths}}(undef, length(aln_asms))
    Threads.@threads for i in 1:length(depths)
        depths[i] =
            load_depths_and_errors(aln_asms[i], paths[i].t_depth, paths[i].a_depth, config)
    end
    segmentorder = map(i -> order_alnasms(i...), zip(aln_asms, depths))

    # Re-order the vectors based on segment and order.
    for i in eachindex(aln_asms, depths, segmentorder)
        ord = sortperm([
            (a.reference.segment, o) for (a, o) in zip(aln_asms[i], segmentorder[i])
        ])
        aln_asms[i] = aln_asms[i][ord]
        depths[i] = depths[i][ord]
        segmentorder[i] = segmentorder[i][ord]
    end

    depth_getters = (i -> i.template_depths, i -> i.assembly_depths)

    # Write depths
    for (path, getter) in zip(
        (joinpath(depths_dir, "template.tsv.gz"), joinpath(depths_dir, "assembly.tsv.gz")),
        depth_getters,
    )
        open(GzipCompressorStream, path, "w") do io
            println(io, "sample\tsegment\torder\tpos\tdepth")
            for (sample, aln_asm_v, depth_v, order_v) in
                zip(samples, aln_asms, depths, segmentorder)
                for (aln_asm, depth, order) in zip(aln_asm_v, depth_v, order_v)
                    segment = aln_asm.reference.segment
                    for (pos, depth_value) in enumerate(getter(depth))
                        println(io, join((sample, segment, order, pos, depth_value), '\t'))
                    end
                end
            end
        end
    end

    passes = report(
        report_path,
        samples,
        segmentorder,
        aln_asms,
        depths,
        read_stats,
        config,
        similar,
    )

    for (p, alnasmv, pass, sorder) in zip(paths, aln_asms, passes, segmentorder)
        write_sequences(p.seq_dir, p.sample, alnasmv, pass, sorder)
    end

    open(GzipCompressorStream, joinpath(tmp_dir, "internal.jls.gz"), "w") do io
        serialize(io, samples, aln_asms, depths, passes, segmentorder)
    end

    # Make plots
    for (sample, aln_asmv, depth_v) in zip(samples, aln_asms, depths)
        for (path, getter) in zip(
            (
                joinpath(depths_dir, "$(sample)_template.pdf"),
                joinpath(depths_dir, "$(sample)_assembly.pdf"),
            ),
            depth_getters,
        )
            Plots.savefig(
                make_depth_plot([
                    (a.reference.segment, getter(d)) for (a, d) in zip(aln_asmv, depth_v)
                ]),
                path,
            )
        end
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
    config::Config,
    similar::Bool,
)::Vector{Vector{Bool}}
    passed = Vector{Bool}[]
    open(report_path, "w") do io
        # Print report for each segment
        for i in eachindex(samples)
            v = report(
                io,
                samples[i],
                order[i],
                alnasms[i],
                depths[i],
                read_stats[i],
                config,
            )
            push!(passed, v)
        end
        print(io, '\n')

        # Check for duplicate segments
        if !similar
            check_duplicates(io, samples, order, alnasms, passed)
        end
    end
    return passed
end

function report(
    io::IO,
    sample::Sample,
    order::Vector{UInt8},
    alnasms::Vector{AlignedAssembly},
    depths::Vector{Depths},
    read_stats::ReadStats,
    config::Config,
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
    sort!(aux; by=i -> first(i).reference.segment)

    # Sample header
    println(io, sample)
    report(io, read_stats)

    # Primary segments report
    for (i, data) in enumerate(primary)
        segment = Segment(i - 1)
        (buf, pass) = report_segment(data, config)
        if data !== nothing
            alnasm, _ = data
            passes[indexof[alnasm]] = pass
        end
        print(io, '\t', pass ? "     " : "FAIL ", rpad(string(segment) * ":", 4))
        write(io, take!(buf))
    end

    # Secondary segments
    if !isempty(aux)
        println(io, "")
        println(io, "\t[ POSSIBLE SUPERINFECTION ]")
        for i in aux
            alnasm, depth, order = i
            _, pass = report_segment((alnasm, depth), config)
            segment = alnasm.reference.segment
            print(
                io,
                "\t",
                pass ? "     " : "FAIL ",
                rpad(string(segment), 3) * ' ' * string(order) * ':',
            )
            print_segment_header(io, alnasm, depth)
            passes[indexof[alnasm]] = pass
        end
    end

    print(io, '\n')
    return passes
end

function report_segment(
    data::Union{Tuple{AlignedAssembly, Depths}, Nothing},
    config::Config,
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
        (error isa Influenza.ErrorLowDepthBases && error.n > config.max_low_depth_bases) || (error isa Influenza.ErrorLowCoverage && error.coverage < 0.75)
    end
    if !isempty(critical_errors)
        println(buf, "\t\tERROR ", first(critical_errors))
        println(buf, "\t\t      Skipping additional errors...")
        return (buf, false)
    end

    for error in alnasm.errors
        p = pass(error, config)
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
            p = pass(error, config)
            println(buf, "\t\t", (p ? "      " : "ERROR "), protein.variant, ": ", error)
            passed &= (p || !is_important(protein.variant))
        end
    end
    return (buf, passed)
end

function print_segment_header(io::IO, alnasm::AlignedAssembly, depth::Depths)
    print(io, " Identity $(@sprintf "%.3f" (alnasm.identity)),")
    println(
        io,
        " depth $(@sprintf "%.2e" assembly_depth(depth)), coverage $(@sprintf "%.3f" template_coverage(depth))",
    )
end

# fallback: segment errors fails
pass(::Influenza.InfluenzaError, ::Config) = false

# Seems weird, but NA stalk can have 75 bp without any problems
pass(x::Influenza.ErrorIndelTooBig, ::Config) = length(x.indel) < 100
pass(x::Influenza.ErrorEarlyStop, ::Config) = x.observed_naa + 14 > x.expected_naa
function pass(x::Influenza.ErrorInsignificant, config::Config)
    x.n_insignificant < ifelse(config.is_illumina, 5, 25)
end
function pass(x::Influenza.ErrorAmbiguous, config::Config)
    x.n_ambiguous < ifelse(config.is_illumina, 5, 25)
end
pass(x::Influenza.ErrorLateStop, ::Config) = true

# We fail with low depth because cross-contamination from other samples
# can give enough reads to reconstruct a segment.
# If more than 500 bp have lower than 25 depth, we can't trust it.
pass(x::Influenza.ErrorLowDepthBases, config::Config) = x.n < config.max_low_depth_bases

function check_duplicates(
    io::IO,
    samples::Vector{Sample},
    order::Vector{Vector{UInt8}},
    alnasms::Vector{Vector{AlignedAssembly}},
    passed::Vector{Vector{Bool}},
)
    kmers = map(zip(alnasms, passed)) do (av, pv)
        map(zip(av, pv)) do (a, p)
            p ? kmerset(a.assembly.seq) : nothing
        end
    end
    println(io, "Possibly duplicated samples:")
    for i in 1:(length(kmers) - 1), j in (i + 1):length(kmers)
        check_duplicates(
            io,
            samples[i],
            samples[j],
            order[i],
            order[j],
            alnasms[i],
            alnasms[j],
            passed[i],
            passed[j],
            kmers[i],
            kmers[j],
        )
    end
end

function check_duplicates(
    io::IO,
    sample1::Sample,
    sample2::Sample,
    order1::Vector{<:Integer},
    order2::Vector{<:Integer},
    alnasms1::Vector{AlignedAssembly},
    alnasms2::Vector{AlignedAssembly},
    passed1::Vector{Bool},
    passed2::Vector{Bool},
    kmers1::Vector{<:Union{Nothing, KMerSet}},
    kmers2::Vector{<:Union{Nothing, KMerSet}},
)
    # Get a list of Bool for each segment, true if at least one segcopy
    # from both samples were passed.
    seg_passed =
        map(&, segments_passed(passed1, alnasms1), segments_passed(passed2, alnasms2))
    n_possible_dup = sum(seg_passed; init=0)

    # We consider at least 2 duplicate segments as a duplicate sample, so if at most
    # 1 segment passed in both samples, we can quit early
    n_possible_dup < 2 && return nothing
    segments_dup = fill(false, N_SEGMENTS)

    pairs = Vector{Tuple{Segment, Float64, UInt8, UInt8}}()
    for i in 1:length(alnasms1), j in 1:length(alnasms2)
        # Only compare passed segcopies
        (passed1[i] & passed2[j]) || continue
        segment = alnasms1[i].reference.segment

        # Of course only compare segcopies of the same segment
        segment == alnasms2[j].reference.segment || continue
        (seq1, seq2) = alnasms1[i].assembly.seq, alnasms2[j].assembly.seq
        len1, len2 = minmax(length(seq1), length(seq2))

        # If the lengths are more than 25% different, they are surely not
        # identical. We do expect some length difference due to primers and such
        len1 / len2 > 0.75 || continue
        k1 = kmers1[i]::KMerSet
        k2 = kmers2[j]::KMerSet

        # The kmer overlap is an inaccurate, but faster way to check for similarity.
        # we use it here to quickly discard pairs which are clearly too different to
        # avoid having to compute an actual alignment
        overlap(k1, k2) ≥ 0.8 || continue
        id = alnid(seq1, seq2)

        # We consider 99.5% identity or above to be the same sequence
        id ≥ 0.995 || continue
        segments_dup[Integer(segment) + 0x01] = true
        push!(pairs, (segment, id, order1[i], order2[j]))
    end
    n_dup = sum(segments_dup; init=0)
    @assert n_possible_dup ≥ n_dup

    # At least two segments must be the same, and at least half the passed segments
    if n_dup < 2 || n_dup < fld(n_possible_dup, 2)
        return nothing
    end
    sort!(pairs; by=first)
    println(io, sample1, ' ', sample2, ": ", n_dup, '/', n_possible_dup, " segments")
    for (segment, id, ord1, ord2) in pairs
        println(
            io,
            '\t',
            rpad(segment, 3),
            " (",
            ord1,
            ") (",
            ord2,
            ") ",
            (@sprintf "%.3f" id),
        )
    end
    println(io)
end

function segments_passed(
    passed::Vector{Bool},
    alnasms::Vector{AlignedAssembly},
)::NTuple{N_SEGMENTS, Bool}
    v = fill(false, N_SEGMENTS)
    for (a, p) in zip(alnasms, passed)
        v[Integer(a.reference.segment) + 0x01] |= p
    end
    NTuple{N_SEGMENTS}(v)
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
        pairalign(OverlapAlignment(), a, b, DEFAULT_DNA_ALN_MODEL),
    )::BioAlignments.PairwiseAlignment
    # And this is only nothing if the identity is very low, but we already
    # checked with kmer overlap, so this should never happen
    Influenza.alignment_identity(OverlapAlignment(), aln)::Float64
end

function make_depth_plot(v::Vector{<:Tuple{Segment, Vector{<:Unsigned}}})
    plt = Plots.plot(; ylabel="Log10 depths", xticks=nothing, ylim=(-0.1, 5))
    for (segment, depth) in v
        ys = log10.(depth)
        xs = range(0.0; stop=1.0, length=length(ys))
        index = Integer(segment) + 1
        Plots.plot!(plt, xs, ys; label=string(segment), legend=:outertopright, color=index)
    end
    return plt
end
