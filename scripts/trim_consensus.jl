# Purpose: Trim primers off consensus seq
# and deduplicate very similar

using FASTX: FASTA
using BioSequences: NucleotideSeq, LongDNASeq, reverse_complement, iscompatible

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
    fuzzylen::Int,
)
    result = seq
    overlaps = [slide!(p, seq, minlength, fuzzylen) for (_, p) in primers]
    if maximum(overlaps) > 0
        i = argmax(overlaps)
        overlap = overlaps[i]
        header = primers[i][1]
        println("Primer $header found in $seqheader with overlap $overlap")
        result = result[(1 + overlap):end]
    end
    return result
end

function remove_primers(
    seqs::Vector{Tuple{String, LongDNASeq}},
    primers::Vector{Tuple{String, LongDNASeq}},
    minlength::Int,
    fuzzylen::Int,
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
)
    primers = load_fna(primerpath)
    # Quick path if no primers to check
    if isempty(primers)
        cp(consensuspath, output)
        return nothing
    end
    consensus = load_fna(consensuspath)
    removed = remove_primers(consensus, primers, minlength, fuzzylen)
    open(FASTA.Writer, output) do writer
        for (header, seq) in removed
            write(writer, FASTA.Record(header, seq))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 5
        error(
            "Usage: julia trim_consensus.jl primers.fna consensus.fna output.fna minlength fuzzylen",
        )
    end
    minlength = parse(Int, ARGS[4])
    minlength < 1 && error("Minlength must be one")
    fuzzylen = parse(Int, ARGS[5])
    fuzzylen >= minlength || error("Fuzzylen cannot be smaller than minlength")

    trim_consensus(ARGS[1], ARGS[2], ARGS[3], minlength, fuzzylen)
end
