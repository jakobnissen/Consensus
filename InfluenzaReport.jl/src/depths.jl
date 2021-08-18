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

function load_depths_and_errors(
    alnasms::Vector{SegmentTuple{Option{AlignedAssembly}}},
    depthspaths::Vector
)::Vector{SegmentTuple{Option{Depths}}}
    depths = load_depths(depthspaths)
    @assert length(alnasms) == length(depthspaths)
    for (alnasm_tup, depth_tup) in zip(alnasms, depths)
        for (m_alnasm, m_depth) in zip(alnasm_tup, depth_tup)
            alnasm = @unwrap_or m_alnasm continue
            depth = @unwrap_or m_depth continue
            add_depths_errors!(alnasm, depth)
        end
    end
    return depths
end

# Single-depthspath method, for approximating assembly-level mapping
function load_depths(
    depthspaths::Vector{String}
)::Vector{SegmentTuple{Option{Depths}}}
    depthspaths |> Map() do depthspath
        matdata = open(GzipDecompressorStream, depthspath) do io
            KMATools.parse_mat(io, depthspath)
        end
        maybe_kmarows_tup = split_mat_segments(matdata, depthspath)
        map(maybe_kmarows_tup) do maybe_kmarows
            and_then(Depths, maybe_kmarows) do kma_rows
                Depths(kma_rows, kma_rows)
            end
        end
    end |> Folds.collect
end

# Same signature as above, except depthspaths are pairs of strings,
# first for template mapping, second for assembly mapping
function load_depths(
    depthspaths::Vector{NamedTuple{(:template, :assembly), Tuple{String, String}}}
)::Vector{SegmentTuple{Option{Depths}}}
    depthspaths |> Map() do depthspair
        temp_tuple = split_mat_segments(open(GzipDecompressorStream, depthspair.template) do io
            KMATools.parse_mat(io, depthspair.template)
        end, depthspair.template)
        assm_tuple = split_mat_segments(open(GzipDecompressorStream, depthspair.assembly) do io
            KMATools.parse_mat(io, depthspair.assembly)
        end, depthspair.assembly)
        ntuple(N_SEGMENTS) do i
            template = @unwrap_or temp_tuple[i] return none(Depths)
            assembly = @unwrap_or assm_tuple[i] return none(Depths)
            some(Depths(template, assembly))
        end
    end |> Folds.collect
end

function split_mat_segments(
    matdata::MatType,
    path::AbstractString
)::SegmentTuple{Option{Vector{KMARowType}}}
    result = fill(none(Vector{KMARowType}), N_SEGMENTS)
    for (header, rows) in matdata
        segment::Segment = let
            p = findlast('_', header)
            s = p === nothing ? nothing : tryparse(Segment, strip(header[p+1:end]))
            s !== nothing ? s : error("Could not parse segment from file $path, header $(header)")
        end
        index = Integer(segment) + 0x01
        is_error(result[index]) || error("Segment $segment present twice in $path")
        result[index] = some(rows)
    end
    return Tuple(result)
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
    dir::AbstractString,
    basenames::Vector{String},
    depthsvec::Vector{SegmentTuple{Option{Depths}}}
)
    isdir(dir) || mkdir(dir)
    @assert length(basenames) == length(depthsvec)
    for (basename, m_depths_tuple) in zip(basenames, depthsvec)
        savefig(
            make_depth_plot(map(m_depths_tuple) do m_depths
                and_then(i -> i.template_depths, Vector{UInt32}, m_depths)
            end),
            joinpath(dir, basename * "_template.pdf")
        )
        if any(m_depths_tuple) do m_depths
            map_or(d -> d.dedicated_assembly, m_depths, false)
        end
            savefig(
                make_depth_plot(map(m_depths_tuple) do m_depths
                    and_then(i -> i.assembly_depths, Vector{UInt32}, m_depths)
                end),
                joinpath(dir, basename * "_assembly.pdf")
            )
        end
    end
end

function make_depth_plot(depths::SegmentTuple{Option{Vector{UInt32}}})
    plt = plot(ylabel="Log10 depths", xticks=nothing, ylim=(-0.1, 5))
    for (index, m_depth) in enumerate(depths)
        depth = @unwrap_or m_depth continue
        segment = Segment(index - 1)
        ys = log10.(depth)
        xs = range(0.0, stop=1.0, length=length(ys))
        plot!(plt, xs, ys, label=string(segment), legend=:outertopright)
    end
    return plt
end