# Explicitly mark sequences as appropriate for phylogenetic analysis.
# E.g. some segments may be failed, but should still be the phylogeny.
# Create a dir "phyloseqs" with sequences.

# The phylogeny pipeline will then preferably read that, or else go in 
# sequences dir if that is not present.

module PhyloSeqs

using Consensus: Consensus, INTERNAL_TYPE
using Influenza: Influenza
using ErrorTypes: unwrap
using FASTX: FASTA

function main(
    internal_path::AbstractString, # path of internal.jls.gz object
    outdir::AbstractString, # directory to create
)
    parent = dirname(rstrip(rstrip(outdir, '/'), '\\'))
    (isempty(parent) || isdir(parent)) || error("Parent directory of output does not exist")
    ispath(outdir) && error("Output path exists: \"$outdir\"")
    internal = Consensus.load_internal(internal_path)
    inds = Consensus.pick_with_preset(i -> i.passed, internal)
    if inds === nothing
        return nothing
    end
    picked = [s for (i, s) in enumerate(internal) if in(i, inds)]
    dump_phyloseqs(outdir, picked)
end

function dump_phyloseqs(outdir::AbstractString, picked::Vector{INTERNAL_TYPE})
    mkdir(outdir)
    by_sample = Dict{String, Vector{INTERNAL_TYPE}}()
    for i in picked
        push!(get!(valtype(by_sample), by_sample, nameof(i.sample)), i)
    end
    for (samplename, internals) in by_sample
        open(FASTA.Writer, joinpath(outdir, samplename * ".fna")) do writer
            for internal in internals
                segment = internal.alnasm.reference.segment
                order = internal.order
                passed = internal.passed ? 'P' : 'F'
                name = "$(samplename)_$(segment)_$(order)_$(passed)"
                seq = internal.alnasm.assembly.seq
                write(writer, FASTA.Record(name, seq))
            end
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 2
        error("Usage: julia internal_jls_path outdir")
    end
    main(ARGS...)
end

end # module
