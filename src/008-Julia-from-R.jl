using DataFrames

# Get a vector of states (xinit) from a dataframe 

function jget_xinit(state::Type{StateXY}, d::DataFrame)
    [StateXY(d.x[i], d.y[i]) for i in 1:nrow(d)]
end

function jget_xinit(state::Type{StateXYZ}, d::DataFrame)
    [StateXYZ(d.x[i], d.y[i], d.z[i]) for i in 1:nrow(d)]
end

function jget_xinit(state::Type{StateXYZD}, d::DataFrame)
    [StateXYZD(d.x[i], d.y[i], d.z[i], d.ang[i]) for i in 1:nrow(d)]
end

# d = DataFrame(x = [1, 2], y = [3, 4])
# jget_xinit(StateXY, d)