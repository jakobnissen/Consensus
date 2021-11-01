"""
    Consensus

This package is the Julia-side (backend) of the flupipe Snakemake pipeline.
It's purpose is to be *application-specific* to the flupipe, NOT to include
generally useful functions.
"""
module Consensus

using FASTX: FASTA
using BioSequences: LongDNASeq, DNA, isgap
using ErrorTypes: Option, some, none, unwrap, unwrap_or, @unwrap_or, and_then, map_or, is_error
using Influenza: Influenza, Segment, Assembly, Reference, AlignedAssembly, Protein
using FluWorkflowTools: Sample, split_segment, try_parseout_suffix
using KMATools: KMATools
using Printf: @sprintf
using Transducers: Map
using Folds: Folds
using Plots: Plots
using CodecZlib: GzipDecompressorStream, GzipCompressorStream
using Serialization: serialize
using JSON3: JSON3

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

const _IMPORTANT = Tuple(Bool[1,1,0,0,1,0,1,1,1,1,1,1,1,0,1,1,1])
@assert length(_IMPORTANT) == length(instances(Protein))
is_important(x::Protein) = @inbounds _IMPORTANT[reinterpret(UInt8, x) + 0x01]

function serialize_alnasms(
    io::IO,
    alnasms::Vector{Vector{AlignedAssembly}},
    passes::Vector{Vector{Bool}},
    order::Vector{Vector{UInt8}},
)
    v = map(zip(alnasms, passes, order)) do (a, p, o)
        collect(zip(a, p, o))
    end
    serialize(io, v)
end

include("readqc.jl")
include("alignedassembly.jl")
include("depths.jl")
include("report.jl")

end # module
