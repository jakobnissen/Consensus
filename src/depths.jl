# Two different notions of depth are relevant: Depth relative to the template mapped against,
# and depth of the final assembly.
# The former is relevant for knowing if the template is closely enough related. This can
# be measured by mapping to the template.
# The latter for knowing if the assembly went well. This can only truly be measured by
# mapping to the assembly - but we don't do that for Nanopore reads, and Medaka doesn't
# give us any good way of measuring depth over the assembly. So we use the template depth
# to estimate the assembly depth.

struct Depths
    # Depth of bases of the template, as given by the initial mapping against it
    template_depths::Vector{UInt32}

    # Depth of bases across the assembly, including any inserted bases and not
    # including deleted bases
    assembly_depths::Vector{UInt32}

    # Mark whether the assembly depths are correctly calculated from the depths
    # over the final assembly, instead of estimated from mapping to the template
    dedicated_assembly::Bool
end

coverage(d::Depths) = count(!iszero, d.template_depths) / length(d.template_depths)
mean_depth(d::Depths) = sum(d.assembly_depths) / length(d.assembly_depths)

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

    return Depths(template_depths, assembly_depths, template_mapping !== assembly_mapping)
end

function load_depths_and_errors(alnasms::Vector{AlignedAssembly}, args...)::Vector{Depths}
    depths = load_depths(alnasms, args...)
    for (alnasm, depth) in zip(alnasms, depths)
        add_depths_errors!(alnasm, depth)
    end
    return depths
end

# Single-depthspath method, for approximating assembly-level mapping
function load_depths(
    alnasms::Vector{AlignedAssembly},
    tdepths::String
)::Vector{Depths}
    map(load_depths_from_mat(alnasms, tdepths)) do rows
        Depths(rows, rows)
    end
end

# Same signature as above, except depthspaths are pairs of strings,
# first for template mapping, second for assembly mapping
function load_depths(
    alnasms::Vector{AlignedAssembly},
    tdepths::String,
    adepths::String
)::Vector{Depths}
    trows = load_depths_from_mat(alnasms, tdepths)
    arows = load_depths_from_mat(alnasms, adepths)
    map(i -> Depths(i...), zip(trows, arows))
end

function get_primary(
    alnasms::Vector{AlignedAssembly},
    depths::Vector{Depths}
)::Vector{Bool}
    primary = Dict{Segment, Depths}()
    max_depth = Dict(s => 0.0 for s in instances(Segment))
    for (d, a) in zip(depths, alnasms)
        segment = a.reference.segment
        mean = mean_depth(d)
        if mean > max_depth[segment]
            max_depth[segment] = mean
            primary[segment] = d
        end
    end
    [d === primary[a.reference.segment] for (d,a) in zip(depths, alnasms)]
end

function load_depths_from_mat(
    alnasms::Vector{AlignedAssembly},
    path::String
)::Vector{Vector{KMARowType}}
    indexof = Dict(
        (a.assembly.name, a.reference.segment) => i
        for (i, a) in enumerate(alnasms)
    )
    data = open(GzipDecompressorStream, path) do io
        KMATools.parse_mat(io, path)
    end
    # There are some segments in kma1.mat.gz which has been removed in
    # earlier states of the pipeline.
    # We set an index of nothing for those, then filter them away
    withindex = map(data) do (header, rows)
        key = split_segment(strip(header))
        (get(indexof, key, nothing), rows)
    end
    filter!(i -> first(i) !== nothing, withindex)
    @assert length(withindex) == length(alnasms)
    map(last, sort!(withindex, by=first))
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
    cov = coverage(depth)
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

function plot_depths(
    template_path::String,
    assembly_path::String,
    v::Vector{Tuple{Segment, Depths}}
)::Nothing
    Plots.savefig(make_depth_plot([(s, d.template_depths) for (s, d) in v]), template_path)
    if isempty(v) || first(v)[2].dedicated_assembly
        Plots.savefig(make_depth_plot([(s, d.assembly_depths) for (s, d) in v]), assembly_path)
    end
    return nothing
end

function make_depth_plot(v::Vector{<:Tuple{Segment, Vector{<:Unsigned}}})
    plt = Plots.plot(ylabel="Log10 depths", xticks=nothing, ylim=(-0.1, 5))
    for (segment, depth) in v
        ys = log10.(depth)
        xs = range(0.0, stop=1.0, length=length(ys))
        index = Integer(segment) + 1
        Plots.plot!(plt, xs, ys, label=string(segment), legend=:outertopright, color=index)
    end
    return plt
end