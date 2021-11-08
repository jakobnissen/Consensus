# In this script, it gathers all good enough hits from the .res file,
# then creates new indexable FASTA files for each samplename based on 
# what segments in .res aligns to

using Influenza: Segment, split_segment
using ErrorTypes
using FASTX: FASTA
using BioSequences: LongDNASeq
import KMATools: KMATools

imap(f) = x -> Iterators.map(f, x)
ifilter(f) = x -> Iterators.filter(f, x)

# This is not optimized for speed, but I seriously doubt it's going to be a problem
function read_res(path::AbstractString)::Vector{String}
    rows = open(io -> KMATools.parse_res(io, path), path)

    bysegment = Dict{Segment, Vector{eltype(rows)}}()
    for row in rows
        _, segment = split_segment(row.template)
        push!(get!(valtype(bysegment), bysegment, segment), row)
    end

    # Criterion by which to select the best segment
    function isworse(x, than)
        # Coverage is most critical: Without it, it's impossible to even make an
        # alignment in the first place. But we can't expect more than around 90%
        # due to deletions, or spurious termini in references etc.
        if min(x.tcov, than.tcov) < 0.9
            return x.tcov < than.tcov
        # Why is identity important? Not sure it needs to be here, but a bad ref
        # could get higher depth due to crossmapping from another segment, e.g. at
        # termini. This factor is small if identity is sufficiently high
        elseif min(x.tid, than.tid) < 0.85
            return x.tid < than.tid
        # If coverage and id is OK, depth is the most important factor, since
        # this is what determines what the majority (consensus) variant is.
        else
            return x.depth < than.depth
        end
    end

    best = Dict(k => first(sort!(v, lt=isworse, rev=true)) for (k, v) in bysegment)

    # For the others (non-best), we have different criteria: We pick all hits
    # with id > 85%, cov > 90% and depth > 0.025 * highest_depth
    # and only if depth > 10
    highest_depth = Dict(k => first(sort(v, by=i -> i.depth, rev=true)) for (k,v) in bysegment)

    result = [b.template for b in values(best)]
    for (segment, v) in bysegment
        vfilt = filter(v[2:end]) do row # exclude best segment at index 1
            row.tcov ≥ 0.90 &&
            row.tid ≥ 0.85 &&
            row.depth > max(10.0, 0.025 * highest_depth[segment].depth)
        end
        append!(result, [i.template for i in vfilt])
    end
    return result
end

function collect_sequences(
        refdir::AbstractString,
        headers::Set{String},
)::Dict{String, FASTA.Record}
    # First get a set of present nums
    records = Dict{String, FASTA.Record}()
    present = Set{String}()
    isempty(headers) && return records
    fastapath = joinpath(refdir, "refs.fna")
    open(FASTA.Reader, fastapath) do reader
        record = FASTA.Record()
        while !eof(reader)
            read!(reader, record)
            header = strip(FASTA.header(record)::String)
            if in(header, headers)
                records[header] = copy(record)
                push!(present, header)
            end
        end
        if headers != present
            miss = setdiff(headers, present)
            error("In FASTA $fastapath following headers missing: $(join(miss, ','))")
        end
    end
    return records
end

function dump_sequences(
    alndir::AbstractString,
    headers::Vector{Tuple{String, Vector{String}}}, # (sample, [headers...])
    records::Dict{String, FASTA.Record}
)
    for (samplename, headers_) in headers
        open(FASTA.Writer, joinpath(alndir, samplename, "cat.fna")) do writer
            for header in headers_
                record = records[header]
                write(writer, record)
            end
        end
    end
end

function main(alndir::AbstractString, refdir::AbstractString)
    headers = readdir(alndir) |> imap() do samplename
        (samplename, read_res(joinpath(alndir, samplename, "initial.res")))
    end |> collect
    all_headers = foldl(headers, init=Set{String}()) do existing, new
        samplename, hs = new
        union!(existing, hs)
    end
    records = collect_sequences(refdir, all_headers)
    dump_sequences(alndir, headers, records)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 2
        error("Usage: julia gather_spa.jl alndir refdir")
    end
    main(ARGS[1], ARGS[2])
end
