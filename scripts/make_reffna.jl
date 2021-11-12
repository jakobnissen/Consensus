# Purpose: Create a FASTA file from an input reference JSON file.

using Influenza: load_references
using FASTX: FASTA

function main(injson::AbstractString, outfna::AbstractString)
    open(FASTA.Writer, outfna) do writer
        for ref in load_references(injson)
            write(writer, FASTA.Record(ref.name, ref.seq))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 2
        error("Usage: julia make_reffna.jl injson outfna")
    end
    main(ARGS[1], ARGS[2])
end
