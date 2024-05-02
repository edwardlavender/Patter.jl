using DataFrames



"""
# Julia from R

A collection of functions that facilitate the translation of inputs from R into Julia.

# Details

* [`julia_get_xinit()`] gets a vector of initial States from a DataFrame;

"""
function julia_get end 

function julia_get_xinit(state::Type{StateXY}, d::DataFrame)
    [StateXY(d.map_value[i], d.x[i], d.y[i]) for i in 1:nrow(d)]
end
@doc (@doc julia_get) julia_get_xinit

function julia_get_xinit(state::Type{StateXYZ}, d::DataFrame)
    [StateXYZ(d.map_value[i], d.x[i], d.y[i], d.z[i]) for i in 1:nrow(d)]
end

function julia_get_xinit(state::Type{StateXYZD}, d::DataFrame)
    [StateXYZD(d.map_value[i], d.x[i], d.y[i], d.z[i], d.angle[i]) for i in 1:nrow(d)]
end

# d = DataFrame(map_value = [0, 1], x = [1, 2], y = [3, 4])
# julia_get_xinit(StateXY, d)