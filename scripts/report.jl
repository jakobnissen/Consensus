using Consensus: Consensus
using CodecZlib: GzipDecompressorStream, GzipCompressorStream
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

function make_depth_tsv(path::String,  data::Vector{Consensus.INTERNAL_TYPE}, depths_getter::Function)
    open(GzipCompressorStream, path, "w") do io
        println(io, "sample\tsegment\torder\tpos\tdepth")
        for (sample, alnasm, depths, passed, order) in data
            segment = alnasm.reference.segment
            for (pos, depth) in enumerate(depths_getter(depths))
                println(io, join((sample, segment, order, pos, depth), '\t'))
            end
        end
    end
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
    if length(ARGS) != 4
        println("Usage: julia report.jl platform selfsimilar out_dir ref_dir")
        exit(1)
    end
    platform, similar_str, outdir, refdir = ARGS
    illumina = if platform == "illumina"
        true
    elseif platform == "nanopore"
        false
    else
        println("Platform must be \"illumina\" or \"nanopore\"")
        exit(1)
    end

    similar = parse(Bool, similar_str)

    reportpath = joinpath(outdir, "report_consensus.txt")
    tmpdir = joinpath(outdir, "tmp")
    alndir = joinpath(tmpdir, "aln")
    consdir = joinpath(outdir, "sequences")
    depthsdir = joinpath(outdir, "depths")
    
    Consensus.snakemake_entrypoint(reportpath, refdir, alndir, consdir, tmpdir, illumina, similar)

    bysample = Dict(
        Sample(name) => Tuple{Segment, Consensus.Depths}[]
        for name in readdir(alndir) if !startswith(name, '.')
    )
    data = open(GzipDecompressorStream, joinpath(tmpdir, "internal.jls.gz")) do io
        deserialize(io)
    end
    for (s, a, d, p, o) in data
        push!(bysample[s], (a.reference.segment, d))
    end
    for (sample, v) in bysample
        tpath = joinpath(depthsdir, nameof(sample) * "_template.pdf")
        apath = joinpath(depthsdir, nameof(sample) * "_assembly.pdf")
        plot_depths(tpath, apath, v)
    end

    # I hate this code pattern but whatever
    make_depth_tsv(joinpath(depthsdir, "template.tsv.gz"), data, x -> x.template_depths)
    make_depth_tsv(joinpath(depthsdir, "assembly.tsv.gz"), data, x -> x.assembly_depths)
end