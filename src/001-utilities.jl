using DataFrames
using Dates

# Check the names of a dataframe included required (req) entries
function check_names(input::DataFrame, req)
    # If req is a single string, convert to array
    req = isa(req, String) ? [req] : req  
    if !all(r -> r ∈ names(input), req)
        cols_missing = filter(r -> !(r ∈ names(input)), req)
        msg = "'input' does not contain all required names. One or more of the following name(s) are missing: '$(join(cols_missing, "', '"))'."
        error(msg)
    end
    return nothing
end

# check_names(DataFrame(x = 1, y = 2), "x")
# check_names(DataFrame(x = 1, y = 2), "z")
# check_names(DataFrame(x = 1, y = 2), ["a", "b"])

# Compute the difference in time between t2 and t1 in seconds
# * `t2` and `t1` may be DateTime or Vector{DateTime}
function diffsecs(t2, t1) 
    Dates.value.(t2 .- t1) ./ 1000
end 

# t1 = DateTime("2024-09-01T12:00:00")
# t2 = DateTime("2024-09-01T12:05:00")
# diffsecs(t2, t1)
# t1 = DateTime("2024-09-01T12:00:00")
# t2 = DateTime("2024-09-01T11:55:00")
# diffsecs(t2, t1)
# t1 = DateTime("2024-09-01T12:00:00")
# t2 = DataFrame(
#     timeline = [DateTime("2024-09-01T12:10:00"), DateTime("2024-09-01T12:15:00")]
# )
# diffsecs(t2.timeline, t1)