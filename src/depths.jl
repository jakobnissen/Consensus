# Two different notions of depth are relevant: Depth relative to the template mapped against,
# and depth of the final assembly.
# The former is relevant for knowing if the template is closely enough related. This can
# be measured by mapping to the template.
# The latter for knowing if the assembly went well. This can only truly be measured by
# mapping to the assembly.

struct Depths
    # Depth of bases of the template, as given by the initial mapping against it
    template_depths::Vector{UInt32}

    # Depth of bases across the assembly, including any inserted bases and not
    # including deleted bases
    assembly_depths::Vector{UInt32}
end

_coverage(v::Vector{<:Unsigned}) = count(!iszero, v) / length(v)
_mean_depth(v::Vector{<:Unsigned}) = sum(v) / length(v)

template_coverage(d::Depths) = _coverage(d.template_depths)
assembly_coverage(d::Depths) = _coverage(d.assembly_depths)
template_depth(d::Depths) = _mean_depth(d.template_depths)
assembly_depth(d::Depths) = _mean_depth(d.assembly_depths)

# This is the type KMATools parses the mat file into
const KMARowType = Tuple{DNA, NTuple{6, UInt32}}
const MatType = Vector{Tuple{String, Vector{KMARowType}}}

function Depths(template_mapping::Vector{KMARowType}, assembly_mapping::Vector{KMARowType})::Depths
    # Template depths: All positions of template mapping, except if 
    # it's a gap in the template
    template_depths = UInt32[]
    for (refbase, rowdepths) in template_mapping
        if !isgap(refbase)
            push!(template_depths, sum(rowdepths))
        end
    end

    # Assembly depth: If there are no reads at all, or the most common base is
    # a deletion
    assembly_depths = UInt32[]
    for (refbase, rowdepths) in assembly_mapping
        if rowdepths[6] != maximum(rowdepths) || iszero(maximum(rowdepths))
            push!(assembly_depths, sum(rowdepths[1:5]))
        end
    end

    return Depths(template_depths, assembly_depths)
end

function load_depths_and_errors(
    alnasms::Vector{AlignedAssembly},
    template_mat_path::AbstractString,
    assembly_mat_path::AbstractString
)::Vector{Depths}
    segment_map = Dict(a.reference.name => a.reference.segment for a in alnasms)
    alnasm_by_refheader = Dict((a.reference.name, a.reference.segment) => a for a in alnasms)
    @assert length(alnasm_by_refheader) == length(alnasms)
    depths = load_depths(
        template_mat_path,
        assembly_mat_path,
        Set(keys(alnasm_by_refheader)),
        segment_map
    )
    for (k, depth) in depths
        alnasm = alnasm_by_refheader[k]
        add_depths_errors!(alnasm, depth)
    end
    return map(last, depths)
end

function load_depths(
    template_mat_path::AbstractString,
    assembly_mat_path::AbstractString,
    headers::Set{Tuple{String, Segment}},
    segment_map::Dict{String, Segment}
)::Vector{Tuple{Tuple{String, Segment}, Depths}}
    tmat = open(GzipDecompressorStream, template_mat_path) do io
        KMATools.parse_mat(io, template_mat_path)
    end
    amat = open(GzipDecompressorStream, assembly_mat_path) do io
        KMATools.parse_mat(io, assembly_mat_path)
    end
    load_depths(template_mat_path, tmat, amat, headers, segment_map)
end

function load_depths(
    tmappath::AbstractString,
    template_mapping::MatType,
    assembly_mapping::MatType,
    headers::Set{Tuple{String, Segment}},
    segment_map::Dict{String, Segment}
)::Vector{Tuple{Tuple{String, Segment}, Depths}}
    # The assembly mapping is done through iterative mapping based on the template
    # mapping, so all segments present in the former must be in the latter, but
    # not necessarily vice versa
    asm_headers = Set(first(i) for i in assembly_mapping)
    temp_headers = Set(first(i) for i in template_mapping)
    if !issubset(asm_headers, temp_headers)
        miss = first(setdiff(asm_headers, temp_headers))
        error("Missing header from template mapping in \"$(tmappath)\": \"$(miss)\"")
    end
    @assert issubset(headers, Set((i, segment_map[i]) for i in asm_headers))
    template_by_header = Dict{String, Vector{KMARowType}}()
    for (header, rows) in template_mapping
        template_by_header[header] = rows
    end
    result = Tuple{Tuple{String, Segment}, Depths}[]
    for (header, rows) in assembly_mapping
        segment = segment_map[header]
        pair = (header, segment)
        pair in headers || continue
        depths = Depths(template_by_header[header], rows)
        push!(result, (pair, depths))
    end
    return result
end

function add_depths_errors!(
    alnasm::AlignedAssembly,
    depth::Depths
)::Nothing
    # Assembly too short
    if length(depth.assembly_depths) < 2*TERMINAL + 1
        push!(alnasm.errors, Influenza.ErrorTooShort(length(depth.assembly_depths)))
    end

    # Too low coverage relative to template
    cov = template_coverage(depth)
    if cov < 0.9
        push!(alnasm.errors, Influenza.ErrorLowCoverage(cov))
    end

    # Too many low coverage positions in assembly
    n_lowdepth = count(i -> i < 25, depth.assembly_depths)
    if !iszero(n_lowdepth)
        push!(alnasm.errors, Influenza.ErrorLowDepthBases(n_lowdepth))
    end
    return nothing
end

"Each instance of a segment is ordered from 1 (primary) to N,
in order to be able to refer to e.g. PB1 # 2 unambiguously."
function order_alnasms(
    alnasms::Vector{AlignedAssembly},
    depths::Vector{Depths}
)::Vector{UInt8}
    bysegment = Dict{Segment, Vector{Tuple{AlignedAssembly, Depths}}}()
    for i in eachindex(alnasms, depths)
        seg = alnasms[i].reference.segment
        push!(get!(valtype(bysegment), bysegment, seg), (alnasms[i], depths[i]))
    end
    result = zeros(UInt8, length(alnasms))
    indexof = Dict(a => i for (i, a) in enumerate(alnasms))
    for (_, v) in bysegment
        sort!(v, rev=true, by=i->assembly_depth(last(i)))
        for (i, (alnasm, _)) in enumerate(v)
            result[indexof[alnasm]] = i
        end
    end
    return result
end