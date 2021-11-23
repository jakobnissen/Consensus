# If this type changes, also change some of the scripts in scripts/, since they
# read and write this serialization. Just grep for INTERNAL_TYPE
const INTERNAL_TYPE = NamedTuple{
    (:sample, :alnasm, :depths, :passed, :order),
    Tuple{Sample, AlignedAssembly, Depths, Bool, UInt8}
}

function serialize(
    io::IO,
    samples::Vector{Sample},
    alnasms::Vector{Vector{AlignedAssembly}},
    depths::Vector{Vector{Depths}},
    passes::Vector{Vector{Bool}},
    order::Vector{Vector{UInt8}},
)
    v = INTERNAL_TYPE[]
    for i in eachindex(samples, alnasms, depths, passes, order)
        for j in eachindex(alnasms[i], depths[i], passes[i], order[i])
            push!(v, (;
                sample=samples[i],
                alnasm=alnasms[i][j],
                depths=depths[i][j],
                passed=passes[i][j],
                order=order[i][j])
            )
        end
    end
    Serialization.serialize(io, v)
end

function load_internal(path::AbstractString)::Vector{INTERNAL_TYPE}
    isfile(path) || error("File not found: ", path)
    open(Serialization.deserialize, GzipDecompressorStream, path)
end

function pick_with_preset(f::Function, segments::Vector{INTERNAL_TYPE})
    selected = findall(f, segments)
    options = map(repr_segment, segments)
    menu = TerminalMenus.MultiSelectMenu(options, pagesize=40, selected=selected)
    println("Pick segments")
    return TerminalMenus.request(menu)
end

function repr_segment(x::INTERNAL_TYPE)
    symbol = x.passed ? '✓' : '✖'
    seg = rpad(string(x.alnasm.reference.segment), 3)
    str = "$(symbol)|$(seg)|$(x.order)|$(nameof(x.sample))"
    return length(str) > 59 ? first(str, 59) * '…' : str
end
