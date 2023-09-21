const ReadStats = NamedTuple{(:bp, :frac_bp_kept, :mean_read_len, :frac_q20, :frac_q30)}

function check_reads(path::AbstractString)::ReadStats
    obj = open(JSON3.read, path)

    bp::Int = obj[:summary][:after_filtering][:total_bases]
    frac_bp_kept::Float64 = bp / obj[:summary][:before_filtering][:total_bases]
    mean_read_len::Float64 = Float64(obj[:summary][:after_filtering][:read1_mean_length])
    r2_len = get(obj[:summary][:after_filtering], :read2_mean_length, nothing)
    if r2_len !== nothing
        mean_read_len = (mean_read_len + r2_len) / 2
    end
    frac_q20::Float64 = obj[:summary][:before_filtering][:q20_rate]
    frac_q30::Float64 = obj[:summary][:before_filtering][:q30_rate]

    return (; bp, frac_bp_kept, mean_read_len, frac_q20, frac_q30)
end

function report(io::IO, stats::ReadStats)
    prefix, corr = if stats.bp ≥ 1_000_000_000
        ("G", stats.bp / 1_000_000_000)
    elseif stats.bp ≥ 1_000_000
        ("M", stats.bp / 1_000_000)
    elseif stats.bp ≥ 1_000
        ("K", stats.bp / 1_000)
    else
        ("", Float64(stats.bp))
    end
    print(io, @sprintf("%.1f", corr), ' ', prefix, "bp, ")
    print(io, @sprintf("%.0f", stats.frac_bp_kept * 100), "% bp kept, ")
    print(io, @sprintf("%.0f", stats.mean_read_len), " bp mean, ")
    print(io, @sprintf("%.0f", stats.frac_q20 * 100), "% Q20, ")
    print(io, @sprintf("%.0f", stats.frac_q30 * 100), "% Q30")
    print(io, '\n')
end
