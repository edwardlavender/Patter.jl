import Base.Threads: @threads

export State, StateXY, StateXYZ, StateCXY, StateCXYZ

"""
    State

`State` is an Abstract Type that defines the animal's state at a given time step. 

# Built-in sub-types

The following sub-types are built-in:

-   `StateXY(map_value, x, y)`: Used for two dimensional (x, y) states;
-   `StateXYZ(map_value, x, y, z)`: Used for three-dimensional (x, y, z) states;
-   `StateCXY(map_value, x, y, heading)`: Used for two-dimensional (x, y) correlated states;
-   `StateCXYZ(map_value, x, y, z, heading)`: Used for three-dimensional (x, y, z) correlated states;

These contain the following fields: 

-   `map_value`: The value of the movement map at coordinates (x, y), required for all `State`s (see [`ModelMove`](@ref));
-   `x`, `y`:  Floats that define the animal's x and y coordinates, required for all `State`s;
-   `z`: A float that defines the animal's z (depth) coordinate, required for all `State`s with a depth component;
-   `heading`: A float that defines the heading, required by `StateCXYZ`;

# Custom sub-types

To define a custom sub-type, such as `StateXYZ`*, simply define a `struct` that is a sub-type of `Patter.State`:

```
struct StateXYZ <: Patter.State
    # Map value
    map_value::Float64
    # Coordinates
    x::Float64
    y::Float64
    z::Float64
end
```

* Note that `StateXYZ` is now implemented.

New states should obey the following requirements:
-   All states must include `map_value`, `x` and `y` fields;
-   For states with a depth dimension, the depth field must be named `z` (for [`Patter.simulate_move()`](@ref));
-   For `R` users, all fields must be of type `Float64` for `Patter.r_get_states()` to parse state vectors;

To use a new `State` sub-type in the simulation of animal movements (via [`simulate_path_walk()`](@ref)) and particle-filtering algorithms, the following steps are also necessary:
-   Define a corresponding [`ModelMove`](@ref) sub-type;
-   (optional) Define `Patter.map_init()` and `Patter.states_init()` methods for [`simulate_states_init()`](@ref) to simulate initial states;
-   Define a [`Patter.simulate_step()`](@ref) method (for [`Patter.simulate_move()`](@ref)) to update the state using a [`ModelMove`](@ref) instance (in [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref));
-   Define a [`Patter.logpdf_step()`](@ref) method (for [`Patter.logpdf_move()`](@ref)) to evaluate the probability density of movement from one state to another (in [`particle_smoother_two_filter()`](@ref));

"""
abstract type State end 

# State sub-types:
# * map_value, x, y are essential (see documentation)
# * This is much faster & facilitates the implementation of the depth observation models

struct StateXY <: State
    # Map value
    map_value::Float64
    # Coordinates
    x::Float64
    y::Float64
end 

# State for RW in XYZ
struct StateXYZ <: State
    # Map value
    map_value::Float64   
    # Coordinates 
    x::Float64
    y::Float64
    z::Float64
end 

# State for CRW in XY
struct StateCXY <: State
    # Map value
    map_value::Float64   
    # Coordinates 
    x::Float64
    y::Float64
    # Horizontal direction
    heading::Float64  
end 

# State for CRW in XYZ
struct StateCXYZ <: State
    # Map value
    map_value::Float64   
    # Coordinates 
    x::Float64
    y::Float64
    z::Float64
    # Horizontal direction
    heading::Float64  
end 


# Determine whether or not a state is valid (in water)
function state_is_valid(state::State, zdim::Bool)
    if zdim
        return is_valid(state.map_value, state.z)
    else 
        return is_valid(state.map_value) 
    end 
end 