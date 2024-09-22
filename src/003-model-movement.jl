using LRUCache: LRU

export ModelMove, ModelMoveXY, ModelMoveXYZD
export simulate_states_init


#########################
#########################
#### ModelMove structures

"""
    Movement models 

[`ModelMove`](@ref) is an Abstract Type that groups movement models. 

# Built-in sub-types

[`ModelMove`](@ref) sub-types define the components of different kinds of movement model. The following sub-types are built-in:

-   `ModelMoveXY(map, mobility, dbn_length, dbn_angle)`: A sub-type for two-dimensional (x, y) random walks, based distributions for step lengths (`dbn_length`) and turning angles (`dbn_angle`);
-   `ModelMoveXYZD(map, mobility, dbn_length, dbn_angle_delta, dbn_z_delta)`: A sub-type for four-dimensional (correlated) random walks, based on distributions for step lengths (`dbn_length`), changes in turning angle (`dbn_angle`) and changes in depth (`dbn_z_delta`);

These contain the following fields: 

-   `map`: A field that defines the arena within which movement occurs. The coordinate reference system of the `map` must align with the other components of the movement model, which typically require a Universal Transverse Mercator (planar) projection with coordinates in metres. `map` is required by all movement models;
-   `mobility`: A number that defines the maximum movement distance between two time steps;
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

New [`ModelMove`](@ref) structures should obey the following requirements:
-   The `map` and `mobility` fields are required by all [`ModelMove`](@ref) sub-types; 
-   By default, `map` is assumed to be a `GeoArray` but a shapefile can be used with a custom [`extract()`](@ref) method;

To use a new [`ModelMove`](@ref) sub-type in the simulation of animal movements (via [`simulate_path_walk()`](@ref)) and particle-filtering algorithms, the following steps are also necessary:
-   Define a corresponding [`State`](@ref) sub-type;
-   (optional) Define a [`Patter.simulate_state_init()`](@ref) method for [`simulate_states_init()`](@ref) to simulate initial states;
-   Define a [`Patter.simulate_step()`](@ref) method (for [`Patter.simulate_move()`](@ref)) to update the state using a [`ModelMove`](@ref) instance (in [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref));
-   Define a [`Patter.logpdf_step()`](@ref) method (for [`Patter.logpdf_move()`](@ref)) to evaluate the probability density of movement from one state to another (in [`two_filter_smoother()`](@ref));

"""
abstract type ModelMove end 

# Random walk in X and Y
struct ModelMoveXY{T, U, V, W} <: ModelMove
    # Environment
    map::T
    # Distribution of step lengths
    mobility::U
    dbn_length::V
    # Distribution of turning angles
    dbn_angle::W
end 
@doc (@doc ModelMove) ModelMoveXY

# Random walk in X, Y and Z
# * TO DO
# @doc (@doc ModelMove) ModelMoveXYZ

# Correlated random walk in X, Y and random walk in Z
struct ModelMoveXYZD{T, U, V, W, X} <: ModelMove
    # Environment
    map::T
    # Step length
    mobility::U
    dbn_length::V
    # Change in turning angle
    # * This must be symmetrically centred around zero
    dbn_angle_delta::W
    # Change in depth 
    dbn_z_delta::X
end 
@doc (@doc ModelMove) ModelMoveXYZD


#########################
#########################
#### Initialise states

"""
    simulate_state_init(state_type::Type{<:State}, model_move::ModelMove, xlim, ylim)

Simulate a (tentative) initial [`State`](@ref) for an animal.

# Arguments:
-   `state_type`: An empty [`State`](@ref) sub-type, such as [`StateXY`](@ref), used for method dispatch only;
-   `model_move`: A [`ModelMove`](@ref) instance;
-   `xlim`, `ylim`: Pairs of numbers that define the boundaries of the area within which initial (`x`, `y`) coordinates are sampled;

# Details

An initial [`State`](@ref) defines the initial (`x`, `y`) location of an animal, plus initial values for any other [`State`](@ref) components (such as the animal's initial bearing, for [`StateXYZD`](@ref)). The purpose of [`simulate_state_init()`](@ref) is to simulate an initial instance of a [`State`](@ref). This is wrapped by the exported function [`simulate_states_init()`](@ref), which simulates a `Vector` of initial [`State`](@ref)s iteratively until `n` valid [`State`](@ref)s are simulated (see [`is_valid()`](@ref)). Initial [`State`](@ref)s are required to initialise the simulation of individual movement paths (e.g., in [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref)). 

[`simulate_state_init()`](@ref) is an internal generic function, with methods for the built-in [`State`](@ref) sub-types. For custom [`State`](@ref) sub-types, write a corresponding [`simulate_state_init()`](@ref) function to use the exported [`simulate_states_init()`](@ref) function. 

# Returns

- A [`State`](@ref) instance;

# See also 

* [`State`](@ref) for [`State`](@ref) sub-types;
* [`extract`](@ref) to extract values from `map` at [`State`](@ref) coordinates;
* [`is_valid()`](@ref) to determine whether or not a simulated state is valid;
* [`simulate_state_init()`](@ref) to simulate an initial [`State`](@ref);
* [`simulate_states_init()`](@ref) to simulate a `Vector` of initial [`State`](@ref)s, iteratively until `n` valid [`State`](@ref)s are simulated;
* [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref) for the front-end functions that use initial [`State`](@ref)s to simulate animal movement paths;

"""
function simulate_state_init end 

function simulate_state_init(state_type::Type{StateXY}, model_move::ModelMove, xlim, ylim)
    x = rand() * (xlim[2] - xlim[1]) + xlim[1]
    y = rand() * (ylim[2] - ylim[1]) + ylim[1]
    map_value = extract(model_move.map, x, y)
    StateXY(map_value, x, y)
end

# function simulate_state_init(state_type::Type{StateXYZ})
# * TO DO

function simulate_state_init(state_type::Type{StateXYZD}, model_move::ModelMove, xlim, ylim)
    x         = rand() * (xlim[2] - xlim[1]) + xlim[1]
    y         = rand() * (ylim[2] - ylim[1]) + ylim[1]
    map_value = extract(model_move.map, x, y)
    z         = rand() * map_value
    angle     = rand() * 2 * pi
    StateXYZD(map_value, x, y, z, angle)
end 

"""
    simulate_states_init(; state_type::Type{<:State}, model_move::ModelMove, n::Int, xlim, ylim)

Simulate a `Vector` of initial [`State`](@ref)s for the simulation of movement paths.

# Arguments (keywords)
- `state_type`: An empty [`State`](@ref) sub-type, such as `StateXY`, used for method dispatch only;
- `model_move`: A [`ModelMove`](@ref) instance;
- `n`: An integer that defines number of initial states to simulate;
- `xlim`, `ylim`: (optional) Pairs of numbers that define the boundaries of the area within initial (`x`, `y`) coordinates are sampled;

# Details

An initial [`State`](@ref) defines the initial (`x`, `y`) location of an animal, plus initial values for any other [`State`](@ref) components (such as the animal's initial bearing, for [`StateXYZD`](@ref)). The purpose of [`simulate_states_init()`](@ref) is to simulate a `Vector` of initial [`State`](@ref) instances. This wraps the internal generic function [`simulate_state_init()`](@ref), iteratively simulating initial [`State`](@ref)s until `n` valid [`State`](@ref)s are simulated (see [`is_valid()`](@ref)). For custom [`State`](@ref) sub-types, a corresponding [`simulate_state_init()`](@ref) method is required to use this function. Initial [`State`](@ref)s are required to initialise the simulation of individual movement paths (e.g., in [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref)). 

# Returns

- A `Vector` of `n` [`State`](@ref) instances;

# See also 

* [`State`](@ref) and [`ModelMove`](@ref) for [`State`](@ref) and movement model sub-types;
* [`extract`](@ref) to extract values from `map` at [`State`](@ref) coordinates;
* [`is_valid()`](@ref) to determine whether or not a simulated state is valid;
* [`simulate_state_init()`](@ref) to simulate an initial [`State`](@ref);
* [`simulate_states_init()`](@ref) to simulate a `Vector` of initial [`State`](@ref)s, iteratively until `n` valid [`State`](@ref)s are simulated;
* [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref) for the front-end functions that use initial [`State`](@ref)s to simulate animal movement paths;

"""
function simulate_states_init(; state_type::Type{<:State}, model_move::ModelMove, n::Int, xlim = nothing, ylim = nothing)

    # (optional) Define xlim and ylim
    bb = GeoArrays.bbox(model_move.map)
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
        pinit = simulate_state_init(state_type, model_move, xlim, ylim)
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
    simulate_step(state::State, model_move::ModelMove, t::Int64)

Simulate a (tentative) step from one location ([`State`](@ref)) into a new location ([`State`](@ref)).

# Arguments

- `state`: A [`State`](@ref) instance that defines the animal's previous [`State`](@ref);
- `model_move`: A [`ModelMove`](@ref) instance;
- `t`: An integer that defines the time step;

# Details

[`simulate_step()`](@ref) is an internal generic function that simulates a new value for the animal's [`State`], that is, the animal's location (and other state components). Methods are provided for the built-in [`State`] and movement model ([`ModelMove`](@ref)) sub-types. For custom [`State`](@ref)s or [`ModelMove`](@ref) sub-types, corresponding methods must be provided. Internally, [`simulate_step()`](@ref) is wrapped by [`simulate_move()`](@ref), which implements[ `simulate_step()`](@ref) iteratively until a valid [`State`](@ref) is simulated (see [`is_valid()`](@ref)). 

# Returns

- A [`State`](@ref) instance;

# See also

* [`State`](@ref) and [`ModelMove`](@ref) for [`State`](@ref) and movement model sub-types;
* [`simulate_step()`](@ref) to simulate a new [`State`](@ref);
* [`is_valid()`](@ref) to determine whether or not a simulated state is valid;
* [`simulate_move()`](@ref) to simulate states iteratively until a valid state is found;
* [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref) for the front-end functions that use these routines to simulate animal movement paths;

"""
function simulate_step end

# RW in X and Y
function simulate_step(state::StateXY, model_move::ModelMoveXY, t::Int64)
    length    = rand(model_move.dbn_length)
    angle     = rand(model_move.dbn_angle)
    x         = state.x + length * cos(angle)
    y         = state.y + length * sin(angle)
    map_value = extract(model_move.map, x, y)
    StateXY(map_value, x, y)
end 

# RW in X, Y and z
# * TO DO

# CRW in X, Y, Z 
function simulate_step(state::StateXYZD, model_move::ModelMoveXYZD, t::Int64)
    length = rand(model_move.dbn_length)
    angle  = state.angle + rand(model_move.dbn_angle_delta)
    x      = state.x + length * cos(angle)
    y      = state.y + length * sin(angle)
    z      = state.z + rand(model_move.dbn_z_delta)
    map_value = extract(model_move.map, x, y)
    StateXYZD(map_value, x, y, z, angle)
end 


"""
    simulate_move(state::State, model_move::ModelMove, t::Int64, n_trial::Real)

Simulate movement from one location ([`State`](@ref)) into a new location ([`State`](@ref)).

# Arguments

- `state`: A [`State`](@ref) instance that defines the animal's previous [`State`](@ref);
- `model_move`: A [`ModelMove`](@ref) instance;
- `t`: An integer that defines the time step;
- `n_trial`: A number that defines the number of attempts to simulate a valid [`State`](@ref);

# Details

[`simulate_move()`](@ref) is an internal function that uses a [`simulate_step()`](@ref) method to simulate new state proposals iteratively until a valid [`State`](@ref) is generated or `n_trial` is reached (see [`is_valid()`](@ref)). For custom [`State`](@ref)s or [`ModelMove`](@ref) sub-types, corresponding [`simulate_step()`](@ref) methods must be provided for this function. [`simulate_move()`](@ref) is used to simulate movement paths (e.g., in [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref)).

# Returns

- A `Tuple` that comprises the simulated [`State`](@ref) instance and the (log) weight; i.e., 
    * `([State], 0.0)` if [`State`](@ref) is valid;
    * `([State], -Inf)` if [`State`](@ref) is invalid;

# See also

* [`State`](@ref) and [`ModelMove`](@ref) for [`State`](@ref) and movement model sub-types;
* [`simulate_step()`](@ref) to simulate a new [`State`](@ref);
* [`is_valid()`](@ref) to determine whether or not a simulated state is valid;
* [`simulate_move()`](@ref) to simulate states iteratively until a valid state is found;
* [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref) for the front-end functions that use these routines to simulate animal movement paths;

"""
function simulate_move(state::State, model_move::ModelMove, t::Int64, n_trial::Real = 100_000)
    # Determine dimension of movement model
    # * The depth dimension (if present) must be named `z`
    zdim  = hasfield(typeof(state), :z)
    # Iteratively sample a valid location in water
    trial = 1
    while trial <= n_trial
        # Simulate a new state
        pstate = simulate_step(state, model_move, t)
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
    logpdf_step(state_from::State, state_to::State, model_move::ModelMove, length, angle)
    
Evaluate the (unnormalised) log probability of an (unrestricted) movement step. 

# Arguments

- `state_from`: A [`State`](@ref) instance that defines a [`State`](@ref) from which the animal moved;
- `state_to`: A [`State`](@ref) instance that defines a [`State`](@ref) into which the animal moved;
- `model_move`: A [`ModelMove`](@ref) instance;
- `t`: An integer that defines the time step;
- `length`: A float that defines the step length (i.e., the Euclidean distance between `state_from` (`x`, `y`) and `state_to` (`x`, `y`));
- `angle`: A float that defines the angle (in polar coordinates) between `state_from` (`x`, `y`) and `state_to` (`x`, `y`);

# Details

[`logpdf_step()`](@ref) is an internal generic function that evaluates the (unnormalised) log probability of an (unrestricted) movement step between two [`State`](@ref)(s) (i.e., locations). Methods are provided for the built-in [`State`](@ref) and [`ModelMove`](@ref) sub-types, but need to be provided for custom sub-types. Internally, [`logpdf_step()`](@ref) is wrapped by [`logpdf_move()`](@ref), which evaluates the log probability of movement between two [`State`](@ref)s, accounting for restrictions to movement; that is, [`logpdf_move()`](@ref) evaluates `logpdf_step(state_from, state_to, model_move, length, angle) + log(abs(determinate)) - log(Z)` where `Z` is the normalisation constant. This is required for particle smoothing (see [`two_filter_smoother()`](@ref)).

# Returns

- A number (log probability); 

# See also 

* [`State`](@ref) and [`ModelMove`](@ref) for [`State`](@ref) and movement model sub-types;
* [`simulate_step()`](@ref) and [`simulate_move()`](@ref) to simulate new [`State`](@ref)s;
* [`logpdf_step()`](@ref) and [`logpdf_move()`](@ref) to evaluate the log-probability of movement between two locations;
* [`logpdf_move_normalisation()`](@ref) for estimation of the normalisation constant;
* [`two_filter_smoother()`](@ref) for the front-end function that uses these routines for particle smoothing;

"""
function logpdf_step end 

function logpdf_step(state_from::StateXY, state_to::StateXY, model_move::ModelMoveXY, t::Int64, length::Float64, angle::Float64) 
    logpdf(model_move.dbn_length, length) + logpdf(model_move.dbn_angle, angle)
end 

# function logpdf_step(state_from::StateXYZ, ...)

function logpdf_step(state_from::StateXYZD, state_to::StateXYZD, model_move::ModelMoveXYZD, t::Int64, length::Float64, angle::Float64) 
    # Compute change in depth 
    z_delta = state_to.z - state_from.z
    # Compute change in angle 
    angle_delta = abs_angle_difference(angle, state_from.angle)
    # Sum up logpdfs 
    logpdf(model_move.dbn_length, length) + logpdf(model_move.dbn_angle_delta, angle_delta) + logpdf(model_move.dbn_z_delta, z_delta)
end 


"""
    logpdf_move(state_from::State, state_to::State, state_zdim::Bool, 
                model_move::ModelMove, t::Int, box, 
                n_sim::Int,
                cache::LRU)
    
Evaluate the log probability of a movement step between two [`State`](@ref)s (`state_from` and `state_to`). 

# Arguments

- `state_from`: A [`State`](@ref) instance that defines a [`State`](@ref) from which the animal moved;
- `state_to`: A [`State`](@ref) instance that defines a [`State`](@ref) into which the animal moved;
- `state_zdim`: A `Boolian` that defines whether or not `state_from` and `state_to` contain a `z` (depth) dimension;
- `model_move`: A [`ModelMove`](@ref) instance;
- `t`: An integer that defines the time step;
- `box`: (optional) A `NamedTuple` (`min_x`, `max_x`, `min_y`, `max_y`) that defines a 'mobility box' within which movements between `state_from` and `state_to` are always (theoretically) legal. This can be provided if `state_from` and `state_to` are [`StateXY`](@ref) instances and `model_move.map` does not contain `NA`s;
- `n_sim`: An integer that defines the number of Monte Carlo simulations (used to approximate the normalisation constant);
- `cache`: A Least Recently Used (LRU) Cache;

# Details

[`logpdf_move()`](@ref) is an internal function that evaluates the log probability of a movement step between two [`State`](@ref)(s) (i.e., locations). This function wraps [`logpdf_step()`](@ref), accounting for accounting for restrictions to movement; that is, [`logpdf_move()`](@ref) evaluates `logpdf_step(state_from, state_to, model_move, length, angle) + log(abs(determinate)) - log(Z)` where `Z` is the normalisation constant. If `state_from` and `state_to` are two-dimensional states (i.e., [`StateXY`](@ref) instances) and `model_move.map` does not contain `NaN`s, a 'mobility `box`' can be provided. This is a `NamedTuple` of coordinates that define the region within which movements between two locations are always theoretically legal. In this instance, the normalisation constant is simply `log(1.0)`. Otherwise, a Monte Carlo simulation of `n_sim` iterations is required to approximate the normalisation constant, accounting for invalid movements, which is more expensive (see [`logpdf_move_normalisation()`](@ref)). [`logpdf_move()`](@ref) is used for particle smoothing (see [`two_filter_smoother()`](@ref)).

# Returns

- A number (log probability); 

# See also 

* [`State`](@ref) and [`ModelMove`](@ref) for [`State`](@ref) and movement model sub-types;
* [`simulate_step()`](@ref) and [`simulate_move()`](@ref) to simulate new [`State`](@ref)s;
* [`logpdf_step()`](@ref) and [`logpdf_move()`](@ref) to evaluate the log-probability of movement between two locations;
* [`logpdf_move_normalisation()`](@ref) for estimation of the normalisation constant;
* [`two_filter_smoother()`](@ref) for the front-end function that uses these routines for particle smoothing;

"""
function logpdf_move(state_from::State, state_to::State, state_zdim::Bool, 
                     model_move::ModelMove, t::Int, box, n_sim::Int = 100, 
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
    # When locations exceed model_move.mobility, we set -Inf density for speed
    # * Consider NearestNeighbors.jl to build a kd tree:
    # * - Provide locations
    # * - Get index of particles within mobility
    # * - This may improve speed by eliminating the need to iterative over all particles in the smoother 
    if model_move.mobility > length
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
        Z = logpdf_move_normalisation(state_from, state_zdim, model_move, t, n_sim, cache)
    end 

    #### Evaluate density 
    logpdf_step(state_from, state_to, model_move, t, length, angle) + log_det - log(Z)

end 


"""
    logpdf_move_normalisation(state::State, state_zdim::Bool, 
                              model_move::ModelMove, t::Int, n_sim::Int, 
                              cache::LRU)

Approximate the normalisation constant for the (log) probability density of movement from one [`State`](@ref) (location) into another. 

# Arguments

- `state_from`: A [`State`](@ref) instance that defines a [`State`](@ref) from which the animal moved;
- `state_zdim`: A `Boolian` that defines whether or not `state_from` contains a `z` (depth) dimension;
- `model_move`: A [`ModelMove`](@ref) instance;
- `t`: An integer that defines the time step;
- `n_sim`: An integer that defines the number of Monte Carlo simulations;
- `cache`: A Least Recently Used (LRU) cache;

# Details

This internal function runs a Monte Carlo simulation of `n_sim` iterations to estimate the normalisation constant for the (log) probability of movement from one [`State`](@ref) (`state_from`) into another. A Beta(1, 1) prior is used to correct for simulations that fail to generate valid move from `state_from`. The normalisation constant for a given [`State`](@ref) is stored in a LRU cache. This function is used by [`logpdf_move()`](@ref) to evaluate the (log) probability of movement between two states, which is required for particle smoothing (see [`two_filter_smoother()`](@ref)).

# Returns 

- A number (the normalisation constant); 

# See also

* [`State`](@ref) and [`ModelMove`](@ref) for [`State`](@ref) and movement model sub-types;
* [`simulate_step()`](@ref) and [`simulate_move()`](@ref) to simulate new [`State`](@ref)s;
* [`logpdf_step()`](@ref) and [`logpdf_move()`](@ref) to evaluate the log-probability of movement between two locations;
* [`logpdf_move_normalisation()`](@ref) for estimation of the normalisation constant;
* [`two_filter_smoother()`](@ref) for the front-end function that uses these routines for particle smoothing;

"""
function logpdf_move_normalisation(state::State, state_zdim::Bool, model_move::ModelMove, t::Int, n_sim::Int)

    # Run simulation 
    k = 0.0
    for i in 1:n_sim
        # Simulate an (unrestricted) step into a new location 
        pstate = simulate_step(state, model_move, t)
        # Validate the step
        if state_is_valid(pstate, state_zdim)
            k += 1.0
        end
    end

    # Evaluate posterior mean
    # * This assumes a Beta(1, 1) prior
    (k + 1) / (n_sim + 2)
end

# Cached version 
function logpdf_move_normalisation(state::State, state_zdim::Bool, model_move::ModelMove, t::Int, n_sim::Int, cache::LRU)
    get!(cache, state) do
        logpdf_move_normalisation(state, state_zdim, model_move, t, n_sim)
    end
end