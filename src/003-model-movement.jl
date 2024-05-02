#########################
#########################
#### ModelMove structures

export ModelMove, ModelMoveXY, ModelMoveXYZD
export simulate_states_init

# ModelMove structures define the parameters of the movement model
# * All ModelMove structures must contain an `env` field that defines the arena within which movement occurs 
# * By default, `env` is assumed to be a GeoArray but a shapefile can be used with a custom extract() method (see spatial.jl)
# * Users can use a provided structure or write their own
# * For custom ModelMoves, new rstep() methods are required


"""
    Movement models 

`ModelMove` is an abstract type that defines the movement model. 

# Subtypes

-   `MoveModelXY`: A subtype for two-dimensional (x, y) random walks, based distributions for step lengths (`dbn_length`) and turning angles (`dbn_angle`);
-   `MoveModelXYZ`: A subtype for three-dimensional (x, y, z) random walks, based on distributions for step lengths (`dbn_length`), turning angles (`dbn_angle`) and changes in depth (`dbn_z_delta`);
-   `MoveModelXYZD`: A subtype for four-dimensional (correlated) random walks, based on distributions for step lengths (`dbn_length`), changes in turning angle (`dbn_angle`) and changes in depth (`dbn_z_delta`);

# Fields

-   `env`: A field field that defines the arena within which movement occurs. This is required by all movement models;
-   `dbn_length`: The distribution of step lengths;
-   `dbn_angle`: The distribution of turning angles;
-   `dbn_angle_delta`: The distribution of changes in turning angle;
-   `dbn_z_delta`: The distribution of changes in depth;

# Details

-   ModelMove structures define the parameters of the movement model;
-   All ModelMove structures must contain an `env` field. 
-   By default, `env` is assumed to be a GeoArray but a shapefile can be used with a custom [`extract()`](@ref) method;
-   Users can use a provided structure or write their own;
-   For custom ModelMoves, new rstep() methods are required;
"""
abstract type ModelMove end 

# Random walk in X and Y
struct ModelMoveXY{T, U, V} <: ModelMove
    # Environment
    env::T
    # Distribution of step lengths
    dbn_length::U
    # Distribution of turning angles
    dbn_angle::V
end 
@doc (@doc ModelMove) ModelMoveXY

# Random walk in X, Y and Z
# * TO DO
@doc (@doc ModelMove) ModelMoveXYZ

# Correlated random walk in X, Y and random walk in Z
struct ModelMoveXYZD{T, U, V, X} <: ModelMove
    # Environment
    env::T
    # Step length
    dbn_length::U 
    # Change in turning angle
    # * This must be symmetrically centred around zero
    dbn_angle_delta::V
    # Change in depth 
    dbn_z_delta::X
end 
@doc (@doc ModelMove) ModelMoveXYZD


#########################
#########################
#### Initialise states

"""
simulate_state_init(state::State, model::MoveModel, xlim, ylim)

Simulate an initial state.

This function is wrapped by the exported function [`simulate_states_init()`](@ref), which simulates a vector of states.

# Arguments:
-   `state_type`: An empty `State` subtype, such as `StateXY`, used for method dispatch only;
-   `move`: A `MoveModel` instance;
-   `xlim`, `ylim`: Pairs of numbers that define the boundaries of the area within which `x` and `y` state values are sampled;

"""
function simulate_state_init end 

function simulate_state_init(state::Type{StateXY}, model::ModelMove, xlim, ylim)
    x = rand() * (xlim[2] - xlim[1]) + xlim[1]
    y = rand() * (ylim[2] - ylim[1]) + ylim[1]
    StateXY(x, y)
end

# function simulate_state_init(state_type::Type{StateXYZ})
# * TO DO

function simulate_state_init(state::Type{StateXYZD}, model::ModelMove, xlim, ylim)
    x     = rand() * (xlim[2] - xlim[1]) + xlim[1]
    y     = rand() * (ylim[2] - ylim[1]) + ylim[1]
    z     = rand() * extract(model.env, x, y)
    angle = rand() * 2 * pi
    StateXYZD(x, y, z, angle)
end 

"""
# Simulate initial states

Simulate a vector of initial states for the simulation of movement paths and the particle filter.

# Arguments
- `state_type`: An empty `State` subtype, such as `StateXY`, used for method dispatch only;
-   `move`: A `MoveModel` instance;
-   `n`: The number of intial states to simulate;
-   `xlim`, `ylim`: (optional) Pairs of numbers that define the boundaries of the area within which `x` and `y` state values are sampled;

"""
function simulate_states_init(state_type, model::ModelMove, n::Int, xlim = nothing, ylim = nothing)

    # (optional) Define xlim and ylim
    env = model.env
    bb = GeoArrays.bbox(env)
    if isnothing(xlim)
        xlim = (bb.min_x, bb.max_x)
    end
    if isnothing(ylim)
        ylim = (bb.min_y, bb.max_y)
    end
    
    # Sample within limits
    # * We assume that there are at least some valid locations within xlim & ylim
    xinit = state_type[]
    while length(xinit) < n
        pinit = simulate_state_init(state_type, model, xlim, ylim)
        valid = is_valid(model.env, pinit.x, pinit.y)
        if valid
            push!(xinit, pinit)
        end 
    end
    
    xinit

end


#########################
#########################
#### Simulate movement steps

# rstep() methods
# * 
# * The user must provide a new method for new State types

"""
    simulate_step(state::State, model::ModelMove, t::Int64)

Simulate a (tentative) step from one location (`State`) into a new location (`State`).

# Details

-   `simulate_step()` is a generic function that simulates a new value for the animal's `state`;
-   Different methods are dispatched according to the `state` and the movement `model`;
-   New methods must be provided for custom states or movement models;
-   Internally, `simulate_step()` is wrapped by `simulate_move()`, which implements `simulate_step()` iteratively until a valid proposal is generated;
"""

# RW in X and Y
function simulate_step(state::StateXY, model::ModelMoveXY, t::Int64)
    length = rand(model.dbn_length)
    angle  = rand(model.dbn_angle)
    x      = state.x + length * cos(angle)
    y      = state.y + length * sin(angle)
    StateXY(x, y)
end 

# RW in X, Y and z
# * TO DO

# CRW in X, Y, Z 
function simulate_step(state::StateXYZD, model::ModelMoveXYZD, t::Int64)
    length = rand(model.dbn_length)
    angle  = state.angle + rand(model.dbn_angle_delta)
    x      = state.x + length * cos(angle)
    y      = state.y + length * sin(angle)
    z      = state.z + rand(model.dbn_z_delta)
    StateXYZD(x, y, z, angle)
end 


"""
    simulate_move(state::State, model::MoveModel, t::Int64, n_trial::Real)

Simulate movement from one location (`State`) into a new location (`State`).

# Details
-   `simulate_move()` is an internal function that uses an [`simulate_step()`](@ref) method to simulate proposals for a new `state` until a valid proposal is generated;

"""
function simulate_move(state::State, model::ModelMove, t::Int64, n_trial::Real = 100_000)
    # Determine dimension of movement model
    # * The depth dimension (if present) must be named `z`
    zdim  = hasfield(typeof(state), :z)
    # Iteratively sample a valid location in water
    trial = 1
    while trial <= n_trial
        # Simulate a new state
        pstate = simulate_step(state, model, t)
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

