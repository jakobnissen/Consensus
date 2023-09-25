using Consensus: Consensus
using JSON3: JSON3

# The plotting is outside the package because it blows up loading times :(
using Plots: Plots

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 5
        println("Usage: julia report.jl platform selfsimilar out_dir ref_dir config_path")
        exit(1)
    end
    platform, similar_str, outdir, refdir, config_path = ARGS
    is_illumina = if platform == "illumina"
        true
    elseif platform == "nanopore"
        false
    else
        println("Platform must be \"illumina\" or \"nanopore\"")
        exit(1)
    end

    # Read config file
    config = Consensus.Config(Dict(JSON3.read(config_path)), is_illumina)
    similar = parse(Bool, similar_str)

    reportpath = joinpath(outdir, "report_consensus.txt")
    tmpdir = joinpath(outdir, "tmp")
    alndir = joinpath(tmpdir, "aln")
    consdir = joinpath(outdir, "sequences")
    depthsdir = joinpath(outdir, "depths")

    Consensus.snakemake_entrypoint(
        reportpath,
        refdir,
        alndir,
        consdir,
        depthsdir,
        tmpdir,
        similar,
        config,
    )
end
