#########################
#########################
#### ModelMove structures

export ModelMove, ModelMoveXY, ModelMoveXYZD

# ModelMove structures define the parameters of the movement model
# * All ModelMove structures must contain an `env` field that defines the arena within which movement occurs 
# * By default, `env` is assumed to be a GeoArray but a shapefile can be used with a custom extract() method (see spatial.jl)
# * Users can use a provided structure or write their own
# * For custom ModelMoves, new rstep() methods are required

abstract type ModelMove end 

# Random walk in X and Y
struct ModelMoveXY{T, U, V} <: ModelMove
    # Environment
    env::T
    # Distribution of step lengths
    dbn_len::U
    # Distribution of turning angles
    dbn_ang::V
end 

# Random walk in X, Y and Z
# * TO DO

# Correlated random walk in X, Y and random walk in Z
struct ModelMoveXYZD{T, U, V, X} <: ModelMove
    # Environment
    env::T
    # Step length
    dbn_len::U 
    # Change in turning angle
    # * This must be symmetrically centred around zero
    dbn_ang_delta::V
    # Change in depth 
    dbn_z_delta::X
end 


#########################
#########################
#### Initialise states

# rinit() methods generate a candidate initial state for pf_init()

function rinit(state::Vector{StateXY}, model::ModelMove, xlim, ylim)
    x = rand() * (xlim[2] - xlim[1]) + xlim[1]
    y = rand() * (ylim[2] - ylim[1]) + ylim[1]
    StateXY(x, y)
end

# function rinit(state::vector{StateXYZ})
# * TO DO

function rinit(state::Vector{StateXYZD}, model::ModelMove, xlim, ylim)
    x = rand() * (xlim[2] - xlim[1]) + xlim[1]
    y = rand() * (ylim[2] - ylim[1]) + ylim[1]
    z = rand() * extract(model.env, x, y)
    ang = rand() * 2 * pi
    StateXYZD(x, y, z, ang)
end 


#########################
#########################
#### Simulate movement steps

# rstep() methods
# * rstep() is an exported function that simulates a new location for a given state/movement model combination
# * The user must provide a new method for new State types

# RW in X and Y
function rstep(state::StateXY, model::ModelMoveXY, t::Int64)
    len = rand(model.dbn_len)
    ang = rand(model.dbn_ang)
    x = state.x + len * cos(ang)
    y = state.y + len * sin(ang)
    StateXY(x, y)
end 

# RW in X, Y and z
# * TO DO

# CRW in X, Y, Z 
function rstep(state::StateXYZD, model::ModelMoveXYZD, t::Int64)
    len = rand(model.dbn_len)
    ang = state.ang + rand(model.dbn_ang_delta)
    x   = state.x + len * cos(ang)
    y   = state.y + len * sin(ang)
    z   = state.z + rand(model.dbn_z_delta)
    StateXYZD(x, y, z, ang)
end 

# rmove()
# * rmove() is an internal function that uses an rstep() method to simulate a new location iteratively
function rmove(state::State, model::ModelMove, t::Int64, n_trial::Int = 1_000)
    # Determine dimension of movement model
    # * The depth dimension (if present) must be named `z`
    zdim  = hasfield(typeof(state), :z)
    # Iteratively sample a valid location in water
    trial = 1
    while trial <= n_trial
        # Simulate a new state
        pstate = rstep(state, model, t)
        # Identify proposal validity 
        if zdim
            valid = is_valid(model.env, pstate.x, pstate.y, pstate.z)
        else 
            valid = is_valid(model.env, pstate.x, pstate.y)
        end 
        # Return valid proposals and the (log) weight (log(1))
        if valid
            return pstate, 0.0
        end 
        trial += 1
    end 

    # For failed iterations, set (log) weight to log(0)
    state, -Inf

end 

