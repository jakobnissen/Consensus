"""
    InfluenzaReport

This package is the Julia-side (backend) of the flupipe Snakemake pipeline.
It's purpose is to be *application-specific* to the flupipe, NOT to include
generally useful functions.
"""
module InfluenzaReport

using FASTX: FASTA
using BioSequences: LongDNASeq, DNA, isgap
using ErrorTypes: Option, some, none, unwrap, unwrap_or, @unwrap_or, and_then, map_or, is_error
using Influenza: Influenza, Segment, Assembly, Reference, AlignedAssembly, Protein
using KMATools: KMATools
using Printf: @sprintf
using Transducers: Map
using Folds: Folds
using Plots: Plots
using CodecZlib: GzipDecompressorStream
using Serialization: serialize

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

function try_split_segment(s::Union{String, SubString{String}})
    p = findlast(isequal(UInt8('_')), codeunits(s))
    p === nothing && return nothing
    seg = tryparse(Segment, SubString(s, p+1:lastindex(s)))
    seg === nothing && return nothing
    return (SubString(s, 1, prevind(s, p)), seg)
end

function split_segment(source::AbstractString, s::Union{String, SubString{String}})
    pair = try_split_segment(s)
    if pair === nothing
        error("In \"" * source * "\", expected NAME_SEGMENT, got \"", s, '\"')
    else
        pair
    end
end

include("alignedassembly.jl")
include("depths.jl")
include("report.jl")

end # module
