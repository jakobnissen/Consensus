# Make reference JSON from internal.jls.gz

module RefJSON

using Consensus: Consensus, INTERNAL_TYPE
using Influenza: Influenza
using ErrorTypes: unwrap

function main(
    inpath::AbstractString,
    outpath::AbstractString
)
    internal = Consensus.load_internal(inpath)
    inds = Consensus.pick_with_preset(i -> i[4], internal)
    picked = [s for (i, s) in enumerate(internal) if in(i, inds)]
    refs = map(Influenza.Reference, picked)
    Influenza.store_references(outpath, refs)
end

function Influenza.Reference(x::INTERNAL_TYPE)
    sample, alnasm = x[1], x[2]
    proteins = map(alnasm.proteins) do protein
        Influenza.ReferenceProtein(
            protein.variant,
            unwrap(protein.orfs)
        )
    end
    Influenza.Reference(
        "$(nameof(sample))_$(alnasm.reference.segment)",
        alnasm.reference.segment,
        alnasm.assembly.seq,
        proteins
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 2
        error("Usage: julia refjls.jl inpath outpath")
    end
    main(ARGS...)
end

end # module