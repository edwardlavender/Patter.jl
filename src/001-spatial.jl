import GeoArrays

export extract

"""
    extract(env::env::GeoArrays.GeoArray, x::Real, y::Real)

Extract the value of a `GeoArray` (such as a bathymetry grid) at a specififed pair of `x` and `y` coordinates.

# Details

In `Patter.jl`, `NaN` elements are taken to define inhospitable habitats, such as land (see `is_valid()`). 

To extract a value from an environmental layer in another format, such as a shapefile, write a custom `extract()` method.

# Returns

- The value of the `GeoArray` for a coordinate pair within the bounds of `env`;
- NaN (of the same type as `env`'s elements) for a coordinate pair beyond the bounds of `env`;
"""
function extract(env::GeoArrays.GeoArray, x::Real, y::Real)
    cell = GeoArrays.indices(env, (x, y))
    if checkbounds(Bool, env, cell)
        env[cell]
    else
        eltype(env)(NaN)
    end
end


"""
    is_valid(env, x::Real, y::Real)
    is_valid(env, x::Real, y::Real, z::Real)

Determine the validity of a point within a given environment.

For the 2D case (`is_valid(env, x, y)`), this function checks if the z-coordinate obtained from `env` for the point (x, y) is not `NaN`.

For the 3D case (`is_valid(env, x, y, z)`), it checks if the z-coordinate obtained for the point (x, y) is not `NaN`, and also ensures that the provided z-coordinate `z` lies within the valid range, specifically 0 < z ≤ z_env.

# Arguments
- `env`: The environment from which z-coordinates are extracted;
- `x`: The x-coordinate of the point;
- `y`: The y-coordinate of the point;
- `z`: The z-coordinate of the point to be validated (only for the 3D case);

# Returns
- `true` if the conditions for validity are met;
- `false` otherwise;
""" 
function is_valid end

# 2D case
function is_valid(env, x::Real, y::Real)
    z_env = extract(env, x, y)
    if isnan(z_env) 
        false
    else
        true
    end 
end

# 3D case
function is_valid(env, x::Real, y::Real, z::Real)
    z_env = extract(env, x, y)
    if isnan(z_env) 
        false
    elseif 0 < z <= z_env
        true 
    else 
        false
    end 
end 


"""
    distance(x0::Real, y0::Real, x1::Real, y1::Real)

Calculate Euclidean distances between coordinates. 

# Arguments
- `x0`, `y0`: The coordinates of the first point;
- `x1`, `y1`: The coordinates of the second point;

# Returns
- A number that defines the distance between two coordinates.

"""
function distance end

# 2D case
function distance(x0::Real, y0::Real, x1::Real, y1::Real)
    # sqrt((x0 - x1)^2 + (y0 - y1)^2)
    hypot(x0 - x1, y0 - y1)
end