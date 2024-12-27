using Crayons
using DataFrames
using Dates


#########################
#########################
#### Console utilities

# Julia warnings: display
# * This function returns 'warning' messages in red
# * In Patter.jl, we use julia_warning() rather than @warn because the former also forces display on Windows 
function julia_warning(msg::String)
    crayon = Crayon(foreground = :red)
    display(crayon("\nWarning: " * msg))
    nothing
end 

# julia_warning("This is a warning!")


#########################
#########################
#### DataFrame utilities

# DataFrame columns: check columns include required entries 
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


#########################
#########################
#### Time utilities 

# Time differences: compute the time difference (s)
# * `t2` and `t1` may be DateTime or Vector{DateTime}
function diffsecs(t2::Union{Dates.DateTime, Vector{Dates.DateTime}}, 
                  t1::Union{Dates.DateTime, Vector{Dates.DateTime}}) 
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

# Time differences: function call duration (s)
function call_duration(call_start::Dates.DateTime)
    return diffsecs(now(), call_start)
end 

# (Internal) check_timeline_*() functions

function check_timeline_entries(t_sim, t_obs)
    issubset(t_obs, t_sim) || error("There are time stamps in `yobs` not found in `timeline`.")
    nothing
end

function check_timeline_spacing(t_sim)
    all(diff(t_sim) .== first(diff(t_sim))) || error("`timeline` must be a sequence of requally spaced time steps.")
    nothing
end

function check_timeline(t_sim, t_obs)
    check_timeline_entries(t_sim, t_obs)
    check_timeline_spacing(t_sim)
    nothing
end