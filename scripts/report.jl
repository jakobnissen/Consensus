fludirname = dirname(dirname(abspath(@__FILE__)))
include(joinpath(fludirname, "InfluenzaReport", "InfluenzaReport.jl"))

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

    reportpath = joinpath(outdir, "report.txt")
    tmpdir = joinpath(outdir, "tmp")
    alndir = joinpath(tmpdir, "aln")
    consdir = joinpath(outdir, "consensus")
    plotdir = joinpath(outdir, "depths")
    
    if illumina
        InfluenzaReport.illumina_snakemake_entrypoint(reportpath, refdir, alndir, consdir, plotdir, tmpdir)
    else
        InfluenzaReport.nanopore_snakemake_entrypoint(reportpath, refdir, alndir, consdir, plotdir, tmpdir)
    end
end
