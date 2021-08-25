# In this script, it gathers all SPA hits, then creates new indexable
# FASTA files for each basename based on what segments is SPA aligns to

using InfluenzaCore
using ErrorTypes
using FASTX
using BioSequences
import KMATools

imap(f) = x -> Iterators.map(f, x)
ifilter(f) = x -> Iterators.filter(f, x)
const N_SEGMENTS = length(instances(Segment))
const SegmentTuple{T} = NTuple{N_SEGMENTS, T}

"Get the highest query coverage if template coverage > 30%, else highest template coverage.
The idea is that if template coverage is 30%, then it's good enough to build an assembly,
and ultimately, query coverage will dominate the consensus anyway. Furthermore, small
levels of contamination will throw off template coverage, but not query coverage."
function readspa(path::AbstractString)::SegmentTuple{Option{UInt}}
    # criteria by which one row is better than another
    isbetter(a, b) = min(a.tcov, b.tcov) > 0.7 ? (a.qcov > b.qcov) : (a.tcov > b.tcov)

    rows = open(io -> KMATools.parse_spa(io, path), path)
    bysegment = Dict{Segment, Vector{eltype(rows)}}()
    for row in rows
        segment = let
            p = findlast('_', row.template)
            s = p === nothing ? nothing : tryparse(Segment, row.template[p+1:end])
            if s === nothing
                error(
                    "In file ", path, " expected header of format \"[NAME]_[SEGMENT]\" got \"",
                    row.template, '"'
                )
            else
                s
            end
        end
        push!(get!(valtype(bysegment), bysegment, segment), row)
    end
    result = fill(none(UInt), N_SEGMENTS)
    for (segment, rows) in bysegment
        result[Integer(segment) + 0x01] = some(first(sort!(rows, lt=isbetter)).num)
    end
    return SegmentTuple(result)
end

function collect_sequences(
        refdir::AbstractString,
        numbers::Vector{Tuple{String, SegmentTuple{Option{UInt}}}}
)::Dict{UInt, FASTA.Record}
    # First get a set of present nums
    records = Dict{UInt, FASTA.Record}()
    present = Set{UInt}()
    for (basename, stuple) in numbers, (i, mnum) in enumerate(stuple)
        num = @unwrap_or mnum continue
        push!(present, num)
    end
    isempty(present) && return records

    fastapath = joinpath(refdir, "refs.fna")
    open(FASTA.Reader, fastapath) do reader
        seqnum = UInt(0)
        record = FASTA.Record()
        lastnum = maximum(present)
        while !eof(reader)
            read!(reader, record)
            seqnum += UInt(1)
            if in(seqnum, present)
                records[seqnum] = copy(record)
            end
            seqnum == lastnum && break
        end
        if seqnum < lastnum
            error("FASTA $fastapath requested num $lastnum, has only $seqnum records")
        end
    end
    return records
end

function dump_sequences(
    alndir::AbstractString,
    numbers::Vector{Tuple{String, SegmentTuple{Option{UInt}}}},
    records::Dict{UInt, FASTA.Record}
)
    for (basename, stuple) in numbers
        writer = open(FASTA.Writer, joinpath(alndir, basename, "cat.fna"))
        for (i, mnum) in enumerate(stuple)
            num = @unwrap_or mnum continue
            record = records[num]
            header = FASTA.identifier(record)::String
            newrecord = FASTA.Record(header, FASTA.sequence(LongDNASeq, record))
            write(writer, newrecord)
        end
        close(writer)
    end
end

function main(alndir::AbstractString, refdir::AbstractString)
    numbers = readdir(alndir) |> imap() do basename
        (basename, readspa(joinpath(alndir, basename, "sparse.spa")))
    end |> collect
    records = collect_sequences(refdir, numbers)
    dump_sequences(alndir, numbers, records)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 2
        error("Usage: julia gather_spa.jl alndir refdir")
    end
    main(ARGS[1], ARGS[2])
end
