using Consensus: Consensus
using CodecZlib: GzipDecompressorStream
using Serialization: deserialize
using Influenza: Sample, Segment

# The plotting is outside the package because it blows up loading times :(
using Plots: Plots

function plot_depths(
    template_path::String,
    assembly_path::String,
    v::Vector{Tuple{Segment, Consensus.Depths}}
)::Nothing
    Plots.savefig(make_depth_plot([(s, d.template_depths) for (s, d) in v]), template_path)
    Plots.savefig(make_depth_plot([(s, d.assembly_depths) for (s, d) in v]), assembly_path)
    return nothing
end

function make_depth_plot(v::Vector{<:Tuple{Segment, Vector{<:Unsigned}}})
    plt = Plots.plot(ylabel="Log10 depths", xticks=nothing, ylim=(-0.1, 5))
    for (segment, depth) in v
        ys = log10.(depth)
        xs = range(0.0, stop=1.0, length=length(ys))
        index = Integer(segment) + 1
        Plots.plot!(plt, xs, ys, label=string(segment), legend=:outertopright, color=index)
    end
    return plt
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 3
        println("Usage: julia report.jl platform out_dir ref_dir")
        exit(1)
    end
    platform, outdir, refdir = ARGS
    illumina = if platform == "illumina"
        true
    elseif platform == "nanopore"
        false
    else
        println("Platform must be \"illumina\" or \"nanopore\"")
        exit(1)
    end

    reportpath = joinpath(outdir, "report_consensus.txt")
    tmpdir = joinpath(outdir, "tmp")
    alndir = joinpath(tmpdir, "aln")
    consdir = joinpath(outdir, "sequences")
    plotdir = joinpath(outdir, "depths")
    
    Consensus.snakemake_entrypoint(reportpath, refdir, alndir, consdir, tmpdir, illumina)

    data = open(GzipDecompressorStream, joinpath(tmpdir, "internal.jls.gz")) do io
        deserialize(io)
    end

    bysample = Dict{Sample, Vector{Tuple{Segment, Consensus.Depths}}}()
    for (s, a, d, p, o) in data
        v = get!(valtype(bysample), bysample, s)
        push!(v, (a.reference.segment, d))
    end
    for (sample, v) in bysample
        tpath = joinpath(plotdir, nameof(sample) * "_template.pdf")
        apath = joinpath(plotdir, nameof(sample) * "_assembly.pdf")
        plot_depths(tpath, apath, v)
    end
end