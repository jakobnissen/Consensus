struct LowQ30
    total_bases::UInt
    q30_bases::UInt

    function LowQ30(_tot::Integer, _q30::Integer)
        tot, q30 = convert(UInt, _tot), convert(UInt, _q30)
        q30 > tot && error("Cannot have more Q30 reads than total reads: $q30 / $tot")
        new(tot, q30)
    end
end

function Base.print(io::IO, x::LowQ30)
    perc = @sprintf "%.3f" (x.q30_bases / x.total_bases)
    print(io, "Reads: Low fraction of Q30: ", perc)
end

struct LowMeanLength
    read1::Float32
    read2::Float32
end

function Base.print(io::IO, x::LowMeanLength)
    print(io, "Reads: Low fw/rv mean len: ", round(Int, x.read1), '/', round(Int, x.read2))
end

struct LowAfterFilter
    total_reads::UInt32
    filter_reads::UInt32

    function LowAfterFilter(_tot::Integer, _filt::Integer)
        tot, filt = convert(UInt32, _tot), convert(UInt32, _filt)
        filt > tot && error("Cannot have more Q30 reads than total reads: $filt / $tot")
        new(tot, filt)
    end
end

function Base.print(io::IO, x::LowAfterFilter)
    perc = @sprintf "%.3f" (x.filter_reads / x.total_reads)
    print(io, "Reads: Low fraction unfiltered: ", perc)
end

const ReadError = Union{LowQ30, LowMeanLength, LowAfterFilter}

function check_reads(path::AbstractString)::Vector{ReadError}
    obj = open(JSON3.read, path)
    result = ReadError[]
    
    total_reads = obj[:summary][:before_filtering][:total_reads]
    filter_reads = obj[:summary][:after_filtering][:total_reads]
    if filter_reads / total_reads < 0.9
        push!(result, LowAfterFilter(total_reads, filter_reads))
    end

    total_bases = obj[:summary][:before_filtering][:total_bases]
    q30_bases = obj[:summary][:before_filtering][:q30_bases]
    if q30_bases / total_bases < 0.9
        push!(result, LowQ30(total_bases, q30_bases))
    end

    read1_length = obj[:summary][:after_filtering][:read1_mean_length]
    read2_length = obj[:summary][:after_filtering][:read2_mean_length]
    if read1_length < 200 || read2_length < 200
        push!(result, LowMeanLength(read1_length, read2_length))
    end

    return result
end