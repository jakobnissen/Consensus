struct Config
    max_low_depth_bases::Int
    depth_threshold::Int
    is_illumina::Bool
end

function Config(dict::Dict{Symbol}, is_illumina::Bool)
    max_low_depth_bases = Int(get(dict, :max_low_depth_bases, 500))::Int
    if max_low_depth_bases < 0
        error("max_low_depth_bases must be at least 0")
    end

    depth_threshold = Int(get(dict, :depth_threshold, 25))::Int
    if max_low_depth_bases < 0
        error("depth_threshold must be at least 0")
    end

    Config(max_low_depth_bases, depth_threshold, is_illumina)
end
