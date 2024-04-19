#########################
#########################
#### State structures

# State structures define the state parameters
# * All states must include `x` and `y`
# * For >=3D states, the depth dimension must be named `z` (for rmove())

export State, StateXY, StateXYZ, StateXYZD

abstract type State end 

# 2D states (2D random walk)
struct StateXY <: State
    # Coordinates
    x::Float64
    y::Float64
end 

# 3D states (3D random walk)
struct StateXYZ <: State
    # Coordinates 
    x::Float64
    y::Float64
    z::Float64
end 

# 4D states (2D CRW and RW in z)
struct StateXYZD <: State
    # Coordinates 
    x::Float64
    y::Float64
    z::Float64
    # Horizontal direction
    ang::Float64  
end 
