function serialize(
    io::IO,
    samples::Vector{Sample},
    alnasms::Vector{Vector{AlignedAssembly}},
    depths::Vector{Vector{Depths}},
    passes::Vector{Vector{Bool}},
    order::Vector{Vector{UInt8}},
)
    v = map(zip(samples, alnasms, depths, passes, order)) do (s, a, d, p, o)
        (nameof(s), collect(zip(a, d, p, o)))
    end

    # If this type changes, also change some of the scripts in scripts/, since they
    # read and write this serialization.
    segdata = Tuple{AlignedAssembly, Depths, Bool, UInt8}
    v::Vector{Tuple{String, Vector{segdata}}}
    Serialization.serialize(io, v)
end
