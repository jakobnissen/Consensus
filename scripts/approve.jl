# Call this from a Consensus output directory, listing the segments that despite
# failing the automated QC should nonetheless be approved.
# Moves them to `primary.fna` if applicable, and updates the pass variable
# in the internal representation

# Usage: julia approve.jl

module Approve

using Consensus: Consensus, INTERNAL_TYPE
using Influenza
using CodecZlib
using Serialization
using REPL.TerminalMenus

function main()
    internal = Consensus.load_internal("tmp/internal.jls.gz")
    chosen = Consensus.pick_with_preset(i -> i[4], internal)
    change_approval!(internal, chosen)
    open(GzipCompressorStream, "tmp/internal.jls.gz", "w") do io
        serialize(io, internal)
    end
    dump_sequences(internal)
end

"Change the approved Bool value in internal type"
function change_approval!(
    internal::Vector{INTERNAL_TYPE},
    chosen::Set{Int},
)
    for (i, (s, a, d, p, o)) in enumerate(internal)
        internal[i] = (s, a, d, in(i, chosen), o)
    end
end

"Use existing Consensus.write_sequences to replace the seq files."
function dump_sequences(internal::Vector{INTERNAL_TYPE})
    bysample = Dict{Sample, Any}()
    for (s, a, d, p, o) in internal
        if !haskey(bysample, s)
            bysample[s] = (
                Influenza.AlignedAssembly[],
                Bool[],
                UInt8[]
            )
        else
            push!(bysample[s][1], a)
            push!(bysample[s][2], p)
            push!(bysample[s][3], o)
        end
    end
    for (sample, tup) in bysample
        path = joinpath("sequences", Influenza.nameof(sample))
        
        # Remove existing paths
        for filename in readdir(path, join=true)
            if endswith(filename, ".fna") || endswith(filename, ".faa")
                rm(filename)
            end
        end

        Consensus.write_sequences(path, sample, tup...)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 0
        error("Usage: julia approve.jl")
    end
    main()
end

end # module