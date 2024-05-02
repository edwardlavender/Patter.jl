export State, StateXY, StateXYZ, StateXYZD

"""
    State

`State` is an abstract type that defines the animal's state at a given time step. 

# Subtypes

-   `StateXY`: Used for two dimensional (x, y) states;
-   `StateXYZ`: Used for three-dimensional (x, y, z) states;
-   `StateXYZD`: Used for four-dimensional (x, y, z, direction) states;

# Fields

-   `x`, `y`:  The animal's x and y coordinates, required for all `States`;
-   `z`: The animal's z coordinate, required for 3D states;
-   `angle`: The turning angle, required by `StateXYZD`;

# Details

-   All states must include `x` and `y` fields;
-   For >= 3D states, the depth dimension must be named `z` (for [`rmove()`](@ref));
-   For `R` users, all fields must be of type `Float64` for [`rget_state_df()`](@ref) to parse state vectors;

# Examples
```jldoctest
julia> StateXY(0.0, 0.0)
StateXY(0.0, 0.0)
```
"""
abstract type State end 

# 2D states (2D random walk)
struct StateXY <: State
    # Coordinates
    x::Float64
    y::Float64
end 
@doc (@doc State) StateXY

# 3D states (3D random walk)
struct StateXYZ <: State
    # Coordinates 
    x::Float64
    y::Float64
    z::Float64
end 
@doc (@doc State) StateXYZ

# 4D states (2D CRW and RW in z)
struct StateXYZD <: State
    # Coordinates 
    x::Float64
    y::Float64
    z::Float64
    # Horizontal direction
    angle::Float64  
end 
@doc (@doc State) StateXYZD

# TO DO
# * Include map_value as the first element of every State and benchmark