# This script contains functionality to make the `refs.fna` and `refs.jls`
# which is used by the consensus pipeline

using Influenza
using FASTX

function make_refs(basepath::AbstractString, v::Vector{Reference})
    store_references(basepath, v)
    open(FASTA.Writer, basepath * ".fna") do writer
        for ref in v
            name = let
                s = string(ref.segment)
                if endswith(ref.name, "_$s")
                    ref.name
                else
                    ref.name * '_' * s
                end
            end
            write(writer, FASTA.Record(name, ref.seq))
        end
    end
end