using Crayons
using DataFrames
using Dates


#########################
#########################
#### Console utilities

# Julia info: display 
# * This function displays info messages 
# * In Patter.jl, we use julia_warning() rather than @warn because the former also forces display on Windows via JuliaCall
function julia_info(msg::String, verbose::Bool = true)
    if verbose
        crayon = Crayon(foreground = :blue)
        display(crayon("Message: " * msg))
    end
    return nothing
end

# Julia warnings: display
# * This function returns 'warning' messages in red
function julia_warning(msg::String)
    crayon = Crayon(foreground = :red)
    display(crayon("\nWarning: " * msg))
    return nothing
end 

# julia_warning("This is a warning!")


#########################
#########################
#### DataFrame utilities

# DataFrame columns: check columns include required entries 
function check_names(input::DataFrame, req::Union{String, Vector{String}})
    # If req is a single string, convert to array
    req = isa(req, String) ? [req] : req  
    if !all(r -> r ∈ names(input), req)
        cols_missing = filter(r -> !(r ∈ names(input)), req)
        msg = "'input' does not contain all required names. One or more of the following name(s) are missing: '$(join(cols_missing, "', '"))'."
        error(msg)
    end
    return nothing
end


#########################
#########################
#### Time utilities 

# Time differences: compute the time difference (s)
# * `t2` and `t1` may be DateTime or Vector{DateTime}
function diffsecs(t2::Union{Dates.DateTime, Vector{Dates.DateTime}}, 
                  t1::Union{Dates.DateTime, Vector{Dates.DateTime}}) 
    Dates.value.(t2 .- t1) ./ 1000
end 

# Time differences: function call duration (s)
function call_duration(call_start::Dates.DateTime)
    return diffsecs(now(), call_start)
end 

# (Internal) check_timeline_*() functions
# * t_obs may be Vector(Any) from R

function check_timeline_entries(t_sim::Vector{Dates.DateTime}, t_obs::Vector)
    issubset(t_obs, t_sim) || error("There are time stamps in `yobs` not found in `timeline`.")
    nothing
end

function check_timeline_spacing(t_sim::Vector{Dates.DateTime})
    all(diff(t_sim) .== first(diff(t_sim))) || error("`timeline` must be a sequence of regularly spaced time steps.")
    nothing
end

function check_timeline(t_sim::Vector{Dates.DateTime}, t_obs::Vector)
    check_timeline_entries(t_sim, t_obs)
    check_timeline_spacing(t_sim)
    nothing
end

#########################
#########################
#### Batching

# Define a Vector of indices for each chunk
# * This function is inspired by parallel::splitIndices() in R
function split_indices(indices::Vector{Int}, n_chunk::Int)
    div, rem = divrem(length(indices), n_chunk)
    splits = [indices[(i-1)*div + min(i-1, rem) + 1 : i*div + min(i, rem)] for i in 1:n_chunk]
    return splits
end 


#########################
#########################
#### Progress

# For R implementations, on Windows, we use a custom DisplayIO for the progress bar (for JuliaCall)
# On MacOS/Linux, we the default stderr
# This is a work around for https://github.com/JuliaInterop/JuliaCall/issues/232

mutable struct ProgressWindowsDisplayIO <: IO
    buffer::String
end

ProgressWindowsDisplayIO() = ProgressWindowsDisplayIO("")

function Base.write(io::ProgressWindowsDisplayIO, s::Union{String, SubString{String}})
    s_clean = replace(s, "\e[K" => "")
    io.buffer *= s_clean
    if occursin('\r', io.buffer) || occursin('\n', io.buffer)
        for part in split(io.buffer, r"[\r\n]") 
            part = strip(part)
            if !isempty(part)
                display(part)
            end
        end
        io.buffer = ""
    end
    return sizeof(s)
end

function Base.write(io::ProgressWindowsDisplayIO, x::UInt8)
    return Base.write(io, string(Char(x)))
end

function Base.write(io::ProgressWindowsDisplayIO, c::Char)
    return Base.write(io, string(c))
end

Base.flush(::ProgressWindowsDisplayIO) = nothing

# Internal function for controlling progress bar from R
# > Functions e.g., particle_filter() have a `progress` argument 
# > This accepts a tuple of arguments passed to Progress
# > For R implementations, we force the output argument to be ProgressWindowsDisplayIO on Windows
function progress_control(; kwargs...)
    if Sys.iswindows()
        return (output = ProgressWindowsDisplayIO(), kwargs...)
    else
        return (; kwargs...)
    end
end