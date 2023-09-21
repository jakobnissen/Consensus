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

function main(path)
    internal = Consensus.load_internal(path)
    chosen = Consensus.pick_with_preset(i -> i.passed, internal)
    if chosen === nothing
        return nothing
    end
    change_approval!(internal, chosen)
    open(GzipCompressorStream, "tmp/internal.jls.gz", "w") do io
        serialize(io, internal)
    end
    dump_sequences(internal)
end

"Change the approved Bool value in internal type"
function change_approval!(internal::Vector{INTERNAL_TYPE}, chosen::Set{Int})
    for (i, x) in enumerate(internal)
        internal[i] = (; x.sample, x.alnasm, x.depths, passed=in(i, chosen), x.order)
    end
end

"Use existing Consensus.write_sequences to replace the seq files."
function dump_sequences(internal::Vector{INTERNAL_TYPE})
    bysample = Dict{Sample, Any}()
    for x in internal
        if !haskey(bysample, x.sample)
            bysample[x.sample] = (Influenza.AlignedAssembly[], Bool[], UInt8[])
        end
        push!(bysample[x.sample][1], x.alnasm)
        push!(bysample[x.sample][2], x.passed)
        push!(bysample[x.sample][3], x.order)
    end
    for (sample, tup) in bysample
        path = joinpath("sequences", Influenza.nameof(sample))

        # Remove existing paths
        for filename in readdir(path; join=true)
            if endswith(filename, ".fna") || endswith(filename, ".faa")
                rm(filename)
            end
        end

        Consensus.write_sequences(path, sample, tup...)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 1
        error("Usage: julia approve.jl internal")
    end
    main(ARGS[1])
end

end # module
