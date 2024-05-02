#########################
#########################
#### State structures

# State structures define the state parameters
# * All states must include `x` and `y`
# * For >= 3D states, the depth dimension must be named `z` (for rmove())
# * All fields must be Float64 for rget_state_df() & correct parsing to R

export State, StateXY, StateXYZ, StateXYZD

"""
    State

The animal's state.
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
    ang::Float64  
end 
@doc (@doc State) StateXYZD