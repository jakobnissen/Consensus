# Purpose: Create a FASTA file from an input reference JSON file.
# And create a .jls serialized file of [(name, segment)] for each ref

using Influenza: load_references
using FASTX: FASTA
using Serialization: serialize

function main(
    injson::AbstractString,
    outfna::AbstractString,
    outjls::AbstractString
)
    refs = load_references(injson)
    segment_map = [(ref.name, ref.segment) for ref in refs]
    open(io -> serialize(io, segment_map), outjls, "w")

    open(FASTA.Writer, outfna) do writer
        for ref in refs
            write(writer, FASTA.Record(ref.name, ref.seq))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 3
        error("Usage: julia make_reffna.jl injson outfna outjls")
    end
    main(ARGS...)
end
