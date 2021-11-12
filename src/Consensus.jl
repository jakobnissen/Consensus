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
using Influenza: Influenza, Sample, Segment, Assembly, Reference, AlignedAssembly, Protein
using KMATools: KMATools
using Printf: @sprintf
using CodecZlib: GzipDecompressorStream, GzipCompressorStream
using Serialization: Serialization
using JSON3: JSON3
using REPL: TerminalMenus

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

include("readqc.jl")
include("alignedassembly.jl")
include("depths.jl")
include("serialization.jl")
include("report.jl")

end # module
