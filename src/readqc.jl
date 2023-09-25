const ReadStats = @NamedTuple{
    bp::Int,
    frac_bp_kept::Float64,
    mean_read_len::Float64,
    frac_q20::Float64,
    frac_q30::Float64,
}

function check_reads(path::AbstractString)::ReadStats
    obj = open(JSON3.read, path)

    bp = Int(obj[:summary][:after_filtering][:total_bases])::Int
    frac_bp_kept = bp / Int(obj[:summary][:before_filtering][:total_bases])::Int
    mean_read_len = Float64(obj[:summary][:after_filtering][:read1_mean_length])::Float64
    r2_len = get(
        obj[:summary][:after_filtering],
        :read2_mean_length,
        nothing,
    )::Union{Int, Nothing}
    if r2_len !== nothing
        mean_read_len = (mean_read_len + r2_len) / 2
    end
    frac_q20 = Float64(obj[:summary][:before_filtering][:q20_rate])::Float64
    frac_q30 = Float64(obj[:summary][:before_filtering][:q30_rate])::Float64

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
