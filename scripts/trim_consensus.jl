# Purpose: Trim primers off consensus seq
# and deduplicate very similar

using FASTX: FASTA
using BioSequences: NucleotideSeq, LongDNASeq, reverse_complement, iscompatible
using Influenza: Influenza, Segment
using FluWorkflowTools: split_segment
using BioAlignments: BioAlignments

const BA = BioAlignments

function is_possibility(needle::NucleotideSeq, haystack::NucleotideSeq, maxerrs::Int)
    errs = 0
    @inbounds for i in eachindex(needle)
        if !iscompatible(needle[i], haystack[i])
            errs += 1
            errs > maxerrs && return false
        end
    end
    return true
end

function slide!(primer::NucleotideSeq, seq::NucleotideSeq, minlen::Int, fuzzylen::Int)
    length(seq) < length(primer) && return 0
    primer = copy(primer)
    seq = copy(seq)[1:length(primer)]
    @assert length(primer) == length(seq)
    for overlap in length(primer):-1:minlen
        maxerrs = ifelse(overlap < fuzzylen, 0, 1)
        is_possibility(primer, seq, maxerrs) && return overlap
        popfirst!(primer)
        pop!(seq)
    end
    return 0
end

function load_fna(path::String)::Vector{Tuple{String, LongDNASeq}}
    open(FASTA.Reader, path) do reader
        map(reader) do record
            h = FASTA.header(record)
            header = isnothing(h) ? "" : h
            (header, FASTA.sequence(LongDNASeq, record))
        end
    end
end

function remove_primers(
    seqheader::String,
    seq::NucleotideSeq,
    primers::Vector{<:Tuple{String, NucleotideSeq}},
    minlength::Int,
    fuzzylen::Int
)
    result = seq
    overlaps = [slide!(p, seq, minlength, fuzzylen) for (_,p) in primers]
    if maximum(overlaps) > 0
        i = argmax(overlaps)
        overlap = overlaps[i]
        header = primers[i][1]
        println("Primer $header found in $seqheader with overlap $overlap")
        result = result[1+overlap:end]
    end
    return result
end

function remove_primers(
    seqs::Vector{Tuple{String, LongDNASeq}},
    primers::Vector{Tuple{String, LongDNASeq}},
    minlength::Int,
    fuzzylen::Int
)::Vector{Tuple{String, LongDNASeq}}
    isempty(seqs) || isempty(primers) && return copy(seqs)
    result = empty(seqs)
    for (header, seq) in seqs
        seq = remove_primers(header, seq, primers, minlength, fuzzylen) # fw
        seq = remove_primers(header, reverse_complement(seq), primers, minlength, fuzzylen) # rv
        seq = reverse_complement(seq)
        push!(result, (header, seq))
    end
    return result
end

function trim_consensus(
    primerpath::String,
    consensuspath::String,
    output::String,
    minlength::Int,
    fuzzylen::Int,
    min_id::Float64
)
    consensus = load_fna(consensuspath)
    primers = load_fna(primerpath)
    removed = remove_primers(consensus, primers, minlength, fuzzylen)
    dedup = deduplicate(removed, consensuspath, min_id)
    open(FASTA.Writer, output) do writer
        for (header, seq) in dedup
            write(writer, FASTA.Record(header, seq))
        end
    end
end

function deduplicate(
    seqs::Vector{Tuple{String, LongDNASeq}},
    path::String,
    min_id::Float64
)::Vector{Tuple{String, LongDNASeq}}
    bysegment = Dict{Segment, Vector{Tuple{String, LongDNASeq}}}()
    for nameseq in seqs
        name, _ = nameseq
        _, segment = split_segment(name)
        push!(get!(valtype(bysegment), bysegment, segment), nameseq)
    end
    result = valtype(bysegment)()
    for (k, v) in bysegment
        append!(result, deduplicate_v(v, min_id))
    end
    return result
end

function deduplicate_v(
    seqs::Vector{Tuple{String, LongDNASeq}},
    min_id::Float64
)::Vector{Tuple{String, LongDNASeq}}
    n = length(seqs)
    n < 2 && return seqs
    kept = [first(seqs)]
    for nameseqi in @view(seqs[2:end])
        _, seqi = nameseqi
        if all(kept) do (_, seqj)
            id_ = id(seqi, seqj)
            id_ === nothing || id_ < min_id
        end
            push!(kept, nameseqi)
        end
    end
    return kept
end

function id(a::LongDNASeq, b::LongDNASeq)
    aln = BA.pairalign(BA.OverlapAlignment(), a, b, Influenza.DEFAULT_DNA_ALN_MODEL).aln
    aln === nothing && return nothing
    return Influenza.alignment_identity(BA.OverlapAlignment(), aln)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 6
        error("Usage: julia trim_consensus.jl primers.fna consensus.fna output.fna minlength fuzzylen id")
    end
    minlength = parse(Int, ARGS[4])
    minlength < 1 && error("Minlength must be one")
    fuzzylen = parse(Int, ARGS[5])
    fuzzylen >= minlength || error("Fuzzylen cannot be smaller than minlength")
    min_id = parse(Float64, ARGS[6])
    if min_id < 0 || min_id > 1
        error("Minimum ID must be in [0.0:1.0]")
    end

    trim_consensus(ARGS[1], ARGS[2], ARGS[3], minlength, fuzzylen, min_id)
end
