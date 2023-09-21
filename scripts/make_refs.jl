"""
    MakeRefs

Create references in a format used by the Consensus pipeline.
"""
module MakeRefs

using Influenza: Influenza, Reference, Segment
using FASTX

imap(f) = x -> Iterators.map(f, x)
ifilter(f) = x -> Iterators.filter(f, x)

struct ReferenceSet
    byid::Dict{String, Reference}
end

ReferenceSet() = ReferenceSet(Dict{String, Reference}())
ReferenceSet(itr) = reduce(push!, itr; init=ReferenceSet())
Base.iterate(x::ReferenceSet, state...) = iterate(values(x.byid), state...)
Base.length(x::ReferenceSet) = length(x.byid)
Base.eltype(::Type{ReferenceSet}) = Reference
Base.copy(x::ReferenceSet) = ReferenceSet(copy(x.byid))

function Base.push!(x::ReferenceSet, y)
    vy = convert(Reference, y)
    x.byid[vy.name] = vy
    return x
end

function Base.append!(x::ReferenceSet, y)
    for i in y
        push!(x, i)
    end
    return x
end

function load_refsets(paths::Vector{<:AbstractString})
    return mapreduce(Influenza.load_references, append!, paths; init=ReferenceSet())
end

"""
Checks the presence of the cd_hit executable.
"""
check_cd_hit() =
    try
        process = run(`cd-hit-est`; wait=false)
        wait(process)
        return true
    catch
        return false
    end

function deduplicate(set::ReferenceSet, id::AbstractFloat, tmpdir=mktempdir())
    check_cd_hit() ||
        error("Command `cd-hit-est` could not be executed. Is the program in your PATH?")
    bysegment = Dict{Segment, Vector{Reference}}()
    for ref in set
        push!(get!(valtype(bysegment), bysegment, ref.segment), ref)
    end

    results = Vector{Vector{Reference}}(undef, length(bysegment))
    Threads.@threads for (i, dat) in collect(enumerate(values(bysegment)))
        results[i] = cd_hit_deduplicate(dat, id, tmpdir)
    end

    return reduce(append!, results; init=ReferenceSet())
end

"Deduplicate using CD-hit"
function cd_hit_deduplicate(
    data::Vector{Reference},
    id::AbstractFloat,
    tmpdir::String=mktempdir(),
)::Vector{Reference}
    # Create FASTA file
    isdir(tmpdir) || error("Directory not found $tmpdir")
    fasta_path, file_io = mktemp(tmpdir)
    writer = FASTA.Writer(file_io)
    for ref in data
        write(writer, FASTA.Record(ref.name, ref.seq))
    end
    close(writer)

    # Run CD hit
    cd_path = run_cd_hit(fasta_path, id)

    # Read results to vector
    names = open(cd_path) do io
        eachline(io) |>
        ifilter(x -> startswith(x, '>')) |>
        imap(x -> String(strip(x)[2:end])) |>
        collect
    end

    byname = Dict(ref.name => ref for ref in data)
    return [byname[i] for i in names]
end

function run_cd_hit(path::AbstractString, id::AbstractFloat)
    outfile = path * ".cdhit"
    command = `cd-hit-est -i $path -o $outfile -aS 0.95 -n 9 -c $(id) -d 32 -T 2 -M 4000`
    pipe = pipeline(command; stdout="$path.log")
    run(pipe)
    return outfile
end

end # module
