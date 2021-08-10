"""
    InfluenzaReport

This package is the Julia-side (backend) of the flupipe Snakemake pipeline.
It's purpose is to be *application-specific* to the flupipe, NOT to include
generally useful functions.
"""
module InfluenzaReport

using FASTX
using BioSequences
using ErrorTypes
using Influenza
using KMATools
using Printf
using Transducers
using Folds
using Plots
using CodecZlib
using Serialization

using Influenza: Assembly, Reference, AlignedAssembly

const N_SEGMENTS = length(instances(Segment))
const SegmentTuple{T} = NTuple{N_SEGMENTS, T}
const TERMINAL = 25

"Too many indel errors in sequence - alignment probably went wrong"
struct ErrorTooManyIndels <: Influenza.ProteinError
    n::UInt32
end

function Base.print(io::IO, x::ErrorTooManyIndels)
    print(io, "Too many indel errors, found ", x.n)
end

# TODO: Superinfection?

const _IMPORTANT = Tuple(Bool[1,1,0,0,1,0,1,1,1,1,1,1,1,0,1,1,1])
@assert length(_IMPORTANT) == length(instances(Protein))
is_important(x::Protein) = @inbounds _IMPORTANT[reinterpret(UInt8, x) + 0x01]

function serialize_alnasms(
    io::IO,
    alnasms::Vector{SegmentTuple{Option{AlignedAssembly}}},
    passes::Vector{SegmentTuple{Bool}}
)
    pairs = map(zip(alnasms, passes)) do (alnasmtuple, passtuple)
        @inbounds ntuple(N_SEGMENTS) do i
            alnasmtuple[i], passtuple[i]
        end
    end
    serialize(io, pairs)
end


include("alignedassembly.jl")
include("depths.jl")
include("report.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 3
        println("Usage: julia InfluenzaReport.jl platform out_dir ref_dir")
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
    alndir = joinpath(outdir, "tmp/aln")
    consdir = joinpath(outdir, "consensus")
    plotdir = joinpath(outdir, "depths")
    tmpdir = joinpath(outdir, "tmp")
    
    if illumina
        illumina_snakemake_entrypoint(reportpath, refdir, alndir, consdir, plotdir, tmpdir)
    else
        nanopore_snakemake_entrypoint(reportpath, refdir, alndir, consdir, plotdir, tmpdir)
    end
end

end # module
