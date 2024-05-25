using GeoArrays

export extract

"""
    extract(map::GeoArrays.GeoArray, x::Real, y::Real)

Extract the value of a `GeoArray` (such as a bathymetry grid) at a pair of `x` and `y` coordinates.

# Details

In `Patter.jl`, `map` is a `GeoArray` that defines the area within which movements are possible. In our applications, `map` is often a bathymetry raster that defines the depth of the seabed across the study area. `map` is an essential field in individual movement models (see [`ModelMove`](@ref)). The [`extract()`](@ref) function supports the simulation of initial states (via [`simulate_states_init()`](@ref)) and the updating of states (via [`simulate_step()`](@ref)), as required to simulate animal movement paths (e.g., in [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref)). Individual states define, at a minimum, the individual's location (`x`, `y`) and the value of the map at that location (`map_value`), which is extracted by [`extract()`](@ref) (see [`State`](@ref)). The Coordinate Reference Systems for `map`, `x` and `y` must align for this to work (i.e., `map` should use a Universal Transverse Mercator projection with coordinates in metres). [`extract()`](@ref) is exported so that it can be used in new methods (for custom states or movement models) of these functions. The simulation of individual movements (via [`simulate_step()`])(@ref) is implemented iteratively (via [`simulate_move()`](@ref)) until a valid movement is found. `NaN` is taken to define inhospitable habitats, such as land, into which the individual cannot move (see [`is_valid()`](@ref)). It should be possible to use a `map` in another format (such as a shapefile) within these routines, with a custom [`extract()`](@ref) method that returns `NaN` in inhospitable habitats and a numeric constant otherwise. 

# See also

* [`State`](@ref) for [`State`](@ref) sub-types;
* [`extract`](@ref) to extract values from `map` at [`State`](@ref) coordinates;
* [`simulate_step()`](@ref) to simulate a new [`State`](@ref);
* [`is_valid()`](@ref) to determine whether or not a simulated state is valid;
* [`simulate_move()`](@ref) to simulate states iteratively until a valid state is found;
* [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref) for the front-end functions that use these routines to simulate animal movement paths;

# Returns

- The value of the `GeoArray` for a coordinate pair within the bounds of `map`;
- NaN (of the same type as `map`'s elements) for a coordinate pair beyond the bounds of `map`;
"""
function extract(map::GeoArrays.GeoArray, x::Real, y::Real)
    # Get row and column indices from coordinates
    rc = GeoArrays.indices(map, [x, y])
    if checkbounds(Bool, map, rc)
        # Extract value of the GeoArray at the specified row/column
        map[rc[1], rc[2], 1]
    else
        eltype(map)(NaN)
    end
end


"""
    is_valid(map_value::Real)
    is_valid(map_value::Real, z::Real)

Determine the validity of a point on a `map`.

For two-dimensional (x, y) states, `is_valid(map_value)` checks if the `map_value` (at a point (`x`, `y`)) is not `NaN`.

For states with a depth (`z`) component, `is_valid(map_value, z)`, checks that `map_value` is not `NaN` and that the provided z-coordinate (`z`) lies within the valid range, specifically `0 < z ≤ map_value` (i.e., the animal is not below the seabed).

# Arguments
- `map_value`: The value of the `map` in a particular location;
- `z`: (optional) The z-coordinate (depth) of the animal, used for states that contain a depth component only;

# Details

These are internal functions. They are used to validate simulated individual states (see [`State`](@ref)). Individual states (i.e., locations) are simulated via [`simulate_step()`](@ref). [`simulate_move()`](@ref) wraps [`simulate_step()`](@ref), iteratively proposing states until a valid state is found. This is required to simulate animal movement paths (e.g., in [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref)).

# See also

* [`State`](@ref) for [`State`](@ref) sub-types;
* [`extract`](@ref) to extract values from `map` at [`State`](@ref) coordinates;
* [`simulate_step()`](@ref) to simulate a new [`State`](@ref);
* [`is_valid()`](@ref) to determine whether or not a simulated state is valid;
* [`simulate_move()`](@ref) to simulate states iteratively until a valid state is found;
* [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref) for the front-end functions that use these routines to simulate animal movement paths;

# Returns
- `true` if the conditions for validity are met;
- `false` otherwise;
""" 
function is_valid end

# 2D case
function is_valid(map_value::Real)
    !isnan(map_value) 
end

# 3D case
function is_valid(map_value::Real, z::Real)
    !isnan(map_value) && 0.0 < z <= map_value
end 

# """
#     bbox(x::Vector{<:Real})

# Define a 'boundary box' `Tuple` from a `Vector` of four numbers (i.e., `min_x`, `max_x`, `min_y` and `max_y`).

# # Details

# The output of this function is designed to match that of [`GeoArray.bbox()`](@ref). This is used to define a 'mobility box' from `R`, for [`logpdf_move()`](@ref). 

# # Returns

# - A `NamedTuple` (`min_x`, `max_x`, `min_y` and `max_y`);

# """  
function bbox(x::Vector{<:Real})
    length(x) == 4 || error("An extent is defined by four numbers.")
    (min_x = x[1], max_x = x[2], min_y = x[3], max_y = x[4])
end 


# """
#     in_bbox(bb, x, y)

# Determine whether or not a coordinate (`x`, `y`) is within a boundary box. This function is used in [`logpdf_move()`](@ref).
# """
function in_bbox(bb, x, y)
    x >= bb.min_x && x <= bb.max_x && y >= bb.min_y && y <= bb.max_y
end


# """
#     distance(x0::Real, y0::Real, x1::Real, y1::Real)

# Calculate Euclidean distances between coordinates. 

# # Arguments
# - `x0`, `y0`: The coordinates of the first point;
# - `x1`, `y1`: The coordinates of the second point;

# # Returns
# - A number that defines the distance between two coordinates.
# """
function distance end

# 2D case
function distance(x0::Real, y0::Real, x1::Real, y1::Real)
    # sqrt((x0 - x1)^2 + (y0 - y1)^2)
    hypot(x0 - x1, y0 - y1)
end


# """
#     cartesian_to_polar(x, y)

# Convert cartesian to polar coordinates. 
# """
function cartesian_to_polar(x, y)
    # (length = sqrt(x^2 + y^2), angle = atan(y, x))
    (length = hypot(x, y), angle = atan(y, x))
end 


# """
#     abs_angle_difference(a1, a2)

# Compute the smallest absolute rotation between two angles.
# """
function abs_angle_difference(a1, a2)
    delta_angle = mod(abs(a1 - a2), 2π)
    min(delta_angle, 2π - delta_angle)
end