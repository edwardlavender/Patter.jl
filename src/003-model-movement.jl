using LRUCache: LRU

export ModelMove, ModelMoveXY, ModelMoveXYZD
export simulate_states_init


#########################
#########################
#### ModelMove structures

"""
    Movement models 

`ModelMove` is an Abstract Type that groups movement models. 

# Built-in sub-types

`ModelMove` sub-types define the components of different kinds of movement model. The following sub-types are built-in:

-   `ModelMoveXY(map, dbn_length, dbn_angle)`: A sub-type for two-dimensional (x, y) random walks, based distributions for step lengths (`dbn_length`) and turning angles (`dbn_angle`);
-   `ModelMoveXYZD(map, dbn_length, dbn_angle_delta, dbn_z_delta)`: A sub-type for four-dimensional (correlated) random walks, based on distributions for step lengths (`dbn_length`), changes in turning angle (`dbn_angle`) and changes in depth (`dbn_z_delta`);

These contain the following fields: 

-   `map`: A field that defines the arena within which movement occurs. The coordinate reference system of the `map` must align with the other components of the movement model, which typically require a Universal Transverse Mercator (planar) projection with coordinates in metres. `map` is required by all movement models;
-   `dbn_length`: The distribution of step lengths;
-   `dbn_angle`: The distribution of turning angles;
-   `dbn_angle_delta`: The distribution of changes in turning angle;
-   `dbn_z_delta`: The distribution of changes in depth;

# Custom sub-types

To define a custom sub-type, such as `ModelMoveXYZ`, simply define a `struct` that is a sub-type of `Patter.ModelMove`:

```
struct ModelMoveXYZ{T, U, V, W} <: Patter.ModelMove
    # The environment (i.e., map)
    # > This defines the regions within which movements are permitted (i.e., in water)
    map::T
    # Distribution for step lengths
    dbn_length::U
    # Distribution for turning angles
    dbn_angle::V
    # Distribution for changes in depth
    dbn_z_delta::W
  end
```

New `ModelMove` structures should obey the following requirements:
-   The `map` field is required by all `ModelMove` sub-types; 
-   By default, `map` is assumed to be a `GeoArray` but a shapefile can be used with a custom [`extract()`](@ref) method;

To use a new `ModelMove` sub-type in the simulation of animal movements (via [`simulate_path_walk()`](@ref)) and particle-filtering algorithms, the following steps are also necessary:
-   Define a corresponding [`State`](@ref) sub-type;
-   (optional) Define a [`Patter.simulate_state_init()`](@ref) method for [`simulate_states_init()`](@ref) to simulate initial states;
-   Define a [`Patter.simulate_step()`](@ref) method (for [`Patter.simulate_move`](@ref)) to update the state using a [`ModelMove`](@ref) instance (in [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref));
-   Define a [`Patter.logpdf_step()`](@ref) method (for [`Patter.logpdf_move`](@ref)) to evaluate the probability density of movement from one state to another (in [`two_filter_smoother()`](@ref));

"""
abstract type ModelMove end 

# Random walk in X and Y
struct ModelMoveXY{T, U, V} <: ModelMove
    # Environment
    map::T
    # Distribution of step lengths
    dbn_length::U
    # Distribution of turning angles
    dbn_angle::V
end 
@doc (@doc ModelMove) ModelMoveXY

# Random walk in X, Y and Z
# * TO DO
# @doc (@doc ModelMove) ModelMoveXYZ

# Correlated random walk in X, Y and random walk in Z
struct ModelMoveXYZD{T, U, V, X} <: ModelMove
    # Environment
    map::T
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
-   `state_type`: An empty `State` sub-type, such as `StateXY`, used for method dispatch only;
-   `model`: A `MoveModel` instance;
-   `xlim`, `ylim`: Pairs of numbers that define the boundaries of the area within which `x` and `y` state values are sampled;

"""
function simulate_state_init end 

function simulate_state_init(state_type::Type{StateXY}, model::ModelMove, xlim, ylim)
    x = rand() * (xlim[2] - xlim[1]) + xlim[1]
    y = rand() * (ylim[2] - ylim[1]) + ylim[1]
    map_value = extract(model.map, x, y)
    StateXY(map_value, x, y)
end

# function simulate_state_init(state_type::Type{StateXYZ})
# * TO DO

function simulate_state_init(state_type::Type{StateXYZD}, model::ModelMove, xlim, ylim)
    x         = rand() * (xlim[2] - xlim[1]) + xlim[1]
    y         = rand() * (ylim[2] - ylim[1]) + ylim[1]
    map_value = extract(model.map, x, y)
    z         = rand() * map_value
    angle     = rand() * 2 * pi
    StateXYZD(map_value, x, y, z, angle)
end 

"""
# Simulate initial states

Simulate a vector of initial states for the simulation of movement paths and the particle filter.

# Arguments
- `state_type`: An empty `State` sub-type, such as `StateXY`, used for method dispatch only;
-   `move`: A `MoveModel` instance;
-   `n`: The number of intial states to simulate;
-   `xlim`, `ylim`: (optional) Pairs of numbers that define the boundaries of the area within which `x` and `y` state values are sampled;

"""
function simulate_states_init(state_type, model::ModelMove, n::Int, xlim = nothing, ylim = nothing)

    # (optional) Define xlim and ylim
    bb = GeoArrays.bbox(model.map)
    if isnothing(xlim)
        xlim = (bb.min_x, bb.max_x)
    end
    if isnothing(ylim)
        ylim = (bb.min_y, bb.max_y)
    end
    
    # Sample within limits
    # * We assume that there are at least some valid locations within xlim & ylim
    xinit = state_type[]
    zdim  = hasfield(state_type, :z)
    while length(xinit) < n
        pinit = simulate_state_init(state_type, model, xlim, ylim)
        valid = state_is_valid(pinit, zdim)
        if valid
            push!(xinit, pinit)
        end 
    end
    
    xinit

end


#########################
#########################
#### Simulate movement steps

# simulate_step() methods
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
function simulate_step end

# RW in X and Y
function simulate_step(state::StateXY, model::ModelMoveXY, t::Int64)
    length    = rand(model.dbn_length)
    angle     = rand(model.dbn_angle)
    x         = state.x + length * cos(angle)
    y         = state.y + length * sin(angle)
    map_value = extract(model.map, x, y)
    StateXY(map_value, x, y)
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
    map_value = extract(model.map, x, y)
    StateXYZD(map_value, x, y, z, angle)
end 


"""
    simulate_move(state::State, model::MoveModel, t::Int64, n_trial::Real)

Simulate movement from one location (`State`) into a new location (`State`).

# Details
-   `simulate_move()` is an internal function that uses an [`Patter.simulate_step()`](@ref) method to simulate proposals for a new `state` until a valid proposal is generated;

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
        valid = state_is_valid(pstate, zdim)
        # Return valid proposals and the (log) weight (log(1))
        if valid
            return pstate, 0.0
        end 
        trial += 1
    end 

    # For failed iterations, set (log) weight to log(0)
    state, -Inf

end 


#########################
#########################
#### Evaluate densities

"""
# Calculate the unnormalised logpdf of an (unrestricted) movement step

# Details

* logpdf_step() is wrapped by the internal logpdf_move() function
* For new states, a corresponding logpdf_step() method is required
* logpdf_move() accounts for restricted steps, the determinate and the normalisation

"""
function logpdf_step(state_from::StateXY, state_to::StateXY, move::ModelMoveXY, length::Float64, angle::Float64) 
    logpdf(move.dbn_length, length) + logpdf(move.dbn_angle, angle)
end 

# function logpdf_step(state_from::StateXYZ, ...)

function logpdf_step(state_from::StateXYZD, state_to::StateXYZD, move::ModelMoveXYZD, length::Float64, angle::Float64) 
    # Compute change in depth 
    z_delta = state_to.z - state_from.z
    # Compute change in angle 
    angle_delta = abs_angle_difference(angle, state_from.angle)
    # Sum up logpdfs 
    logpdf(move.dbn_length, length) + logpdf(move.dbn_angle_delta, angle_delta) + logpdf(move.dbn_z_delta, z_delta)
end 

"""
# Calculate the logpdf of a restricted movement

An internal function that calculates the logpdf of a movement from one state to another. 
"""
function logpdf_move(state_from::State, state_to::State, state_zdim::Bool, 
                     move::ModelMove, t::Int, box, nMC::Int = 100, 
                     cache::LRU = LRU{eltype(state_from), Float64}(maxsize = 100)) 

    #### Validate state 
    # (This check is cheap)
    if !state_is_valid(state_to, state_zdim)
        return -Inf
    end

    #### Evaluate step lengths and angles 
    x = state_to.x - state_from.x
    y = state_to.y - state_from.y
    length, angle = cartesian_to_polar(x, y)
    # When locations are far apart, we set -Inf density for speed
    if cdf(move.dbn_length, length) > 0.999
        return -Inf
    end
    # Calculate log(abs(determinate))
    log_det = -0.5 * log(x^2 + y^2)
    
    #### Approximate the normalisation constant
    # (A) Use box 
    # * Set normalisation constant to 1.0 if the individual is within a box
    # * `box` can be provided for 2D states & if there are no NaNs in the study area
    # * It defines the box within which a move is always valid 
    # * (i.e., the extent of the study area - mobility)
    # * This code is implemented outside of logpdf_move_normalisation(): caching with the box is MUCH slower
    if !isnothing(box) && in_bbox(box, state_from.x, state_from.y)
        Z = 1.0
    # (B) Run MC simulation 
    # * Run simulation
    # * Caching is MUCH faster when MC iterations are required
    else 
        Z = logpdf_move_normalisation(state_from, state_zdim, move, t, nMC, cache)
    end 

    #### Evaluate density 
    logpdf_step(state_from, state_to, move, length, angle) + log_det - log(Z)

end 


# (Internal) normalisation constants
function logpdf_move_normalisation(state::State, state_zdim::Bool, move::ModelMove, t::Int, nMC::Int)

    # Run simulation 
    k = 0.0
    for i in 1:nMC
        # Simulate an (unrestricted) step into a new location 
        pstate = simulate_step(state, move, t)
        # Validate the step
        if state_is_valid(pstate, state_zdim)
            k += 1.0
        end
    end

    # Evaluate posterior mean
    # * This assumes a Beta(1, 1) prior
    (k + 1) / (nMC + 2)
end

# Cached version 
function logpdf_move_normalisation(state::State, state_zdim::Bool, move::ModelMove, t::Int, nMC::Int, cache::LRU)
    get!(cache, state) do
        logpdf_move_normalisation(state, state_zdim, move, t, nMC)
    end
end