using GeoArrays

export extract

"""
    extract(map::GeoArrays.GeoArray, x::Real, y::Real)

Extract the value of a `GeoArray` (such as a bathymetry grid) at a specififed pair of `x` and `y` coordinates.

# Details

In `Patter.jl`, `NaN` elements are taken to define inhospitable habitats, such as land (see `is_valid()`). 

To extract a value from a map in another format, such as a shapefile, write a custom `extract()` method.

# Returns

- The value of the `GeoArray` for a coordinate pair within the bounds of `map`;
- NaN (of the same type as `map`'s elements) for a coordinate pair beyond the bounds of `map`;
"""
function extract(map::GeoArrays.GeoArray, x::Real, y::Real)
    # Get row and column indices from coordinates
    rc = GeoArrays.indices(map, [x, y])
    if checkbounds(Bool, map, rc)
        # Extract value of the GeoArray at the specified row/column
        map[rc[1], rc[1]]
    else
        eltype(map)(NaN)
    end
end


"""
    is_valid(map_value::Real)
    is_valid(map_valu::Real, z::Real)

Determine the validity of a point on a map.

For 2D (x, y) states, `is_valid(map_value)` checks if the `map_value` at point (x, y) is not `NaN`.

For 3D (x, y, z) states, `is_valid(map_value, z)`, checks that `map_value` at point (x, y) is not `NaN`, and also ensures that the provided z-coordinate (`z`) lies within the valid range, specifically 0 < z ≤ `map_value`.

# Arguments
- `map_value`: The map from which z-coordinates are extracted;
- `x`: The x-coordinate of the point;
- `y`: The y-coordinate of the point;
- `z`: The z-coordinate of the point to be validated (only for the 3D case);

# Returns
- `true` if the conditions for validity are met;
- `false` otherwise;
""" 
function is_valid end

# 2D case
function is_valid(map_value::Real)
    if isnan(map_value) 
        false
    else
        true
    end 
end

# 3D case
function is_valid(map_value::Real, z::Real)
    if isnan(map_value) 
        false
    elseif 0.0 < z <= map_value
        true 
    else 
        false
    end 
end 


"""
Determine whether or not a coordinate (x, y) is within a boundary box.
"""
function in_bbox(bb, x, y)
    x >= bb.min_x && x <= bb.max_x && y >= bb.min_y && y <= bb.max_y
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


"""
# Cartesian to polar coordinates
"""
function cartesian_to_polar(x, y)
    (length = hypot(x, y), angle = atan(y, x))
end 

"""
Compute the smallest absolute rotation between two angles
"""
function abs_angle_difference(a1, a2)
    delta_angle = mod(abs(a1 - a2), 2π)
    min(delta_angle, 2π - delta_angle)
end