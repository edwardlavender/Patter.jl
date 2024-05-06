export State, StateXY, StateXYZ, StateXYZD

"""
    State

`State` is an abstract type that defines the animal's state at a given time step. 

# Subtypes

-   `StateXY`: Used for two dimensional (x, y) states;
-   `StateXYZ`: Used for three-dimensional (x, y, z) states;
-   `StateXYZD`: Used for four-dimensional (x, y, z, direction) states;

# Fields

-   `map_value`: The value of the map at coordinates (x, y), required for All `State`s;
-   `x`, `y`:  The animal's x and y coordinates, required for all `State`s;
-   `z`: The animal's z coordinate, required for 3D `State`s;
-   `angle`: The turning angle, required by `StateXYZD`;

# Details

-   All states must include `x`, `y` and `map_value` fields;
-   For >= 3D states, the depth dimension must be named `z` (for [`Patter.simulate_move()`](@ref));
-   For `R` users, all fields must be of type `Float64` for [`Patter.r_get_states()`](@ref) to parse state vectors;
```
"""
abstract type State end 

# State subtypes:
# * map_value, x, y are essential (see documentation)
# * This is much faster & facilitates the implementation of the depth observation models

# 2D states (2D random walk)
struct StateXY <: State
    # Map value
    map_value::Float64
    # Coordinates
    x::Float64
    y::Float64
end 
@doc (@doc State) StateXY

# 3D states (3D random walk)
struct StateXYZ <: State
    # Map value
    map_value::Float64   
    # Coordinates 
    x::Float64
    y::Float64
    z::Float64
end 
@doc (@doc State) StateXYZ

# 4D states (2D CRW and RW in z)
struct StateXYZD <: State
    # Map value
    map_value::Float64   
    # Coordinates 
    x::Float64
    y::Float64
    z::Float64
    # Horizontal direction
    angle::Float64  
end 
@doc (@doc State) StateXYZD


"""
    state_is_valid(state::State, zdim::Bool)

Determine whether or not a `state` is valid.

See also [`State`](@ref), [`is_valid()`](@ref), 
"""
function state_is_valid(state::State, zdim::Bool)
    if zdim
        return is_valid(state.map_value, state.z)
    else 
        return is_valid(state.map_value) 
    end 
end 

