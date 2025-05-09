export ModelMove, ModelMoveXY, ModelMoveXYZ, ModelMoveCXY, ModelMoveCXYZ


#########################
#########################
#### ModelMove structures

"""
    Movement models 

[`ModelMove`](@ref) is an Abstract Type that groups movement models. 

# Built-in sub-types

[`ModelMove`](@ref) sub-types define the components of different kinds of movement model. The following sub-types are built-in:

-   Random walks (RWs):
        - `ModelMoveXY(map, mobility, dbn_length, dbn_heading)`: A sub-type for two-dimensional (x, y) random walks, based distributions for step lengths (`dbn_length`) and headings (`dbn_heading`);
        -  `ModelMoveXYZ(map, mobility, dbn_length, dbn_heading, dbn_z)`: A subtype for a RW in X, Y and Z, based on distributions for step lengths (`dbn_length`), headings (`dbn_heading`) and depth (`dbn_z`);

-   Correlated random walks (CRWs):
        -   `ModelMoveCXY(map, mobility, dbn_length, dbn_heading_delta)`: A sub-type for a CRW in X and Y, based on distributions for step lengths (`dbn_length`) and turning angles (`dbn_heading_delta`) ;
        -   `ModelMoveCXYZ(map, mobility, dbn_length, dbn_heading_delta, dbn_z_delta)`: A sub-type for a CRW in X, Y and Z, based on distributions for step lengths (`dbn_length`), turning angles (`dbn_heading_delta`) and changes in depth (`dbn_z_delta`);

These contain the following fields: 

-   `map`: A field that defines the arena within which movement occurs. The coordinate reference system of the `map` must align with the other components of the movement model, which typically require a Universal Transverse Mercator (planar) projection with coordinates in metres. `map` is required by all movement models;
-   `mobility`: A number that defines the maximum movement distance between two time steps;
-   `dbn_length`: The distribution of step lengths;
-   `dbn_heading`: The distribution of headings;
-   `dbn_heading_delta`: The distribution of turning angles;
-   `dbn_z`: The distribution of depths;
-   `dbn_z_delta`: The distribution of changes in depth;

# Custom sub-types

To define a custom sub-type, such as `ModelMoveXYZ`*, simply define a `struct` that is a sub-type of `Patter.ModelMove`:

```
struct ModelMoveXYZ{T, U, V, W, X} <: Patter.ModelMove
    # The environment (i.e., map)
    # > This defines the regions within which movements are permitted (i.e., in water)
    map::T
    # Distribution for step lengths
    mobility::U
    dbn_length::V
    # Distribution for headings
    dbn_heading::W
    # Distribution for changes in depth
    dbn_z::X
  end
```

* Note this structure is now implemented. 

New [`ModelMove`](@ref) structures should obey the following requirements:
-   The `map` and `mobility` fields are required by all [`ModelMove`](@ref) sub-types; 
-   By default, `map` is assumed to be a `GeoArray` but a shapefile can be used with a custom [`extract()`](@ref) method;

To use a new [`ModelMove`](@ref) sub-type in the simulation of animal movements (via [`simulate_path_walk()`](@ref)) and particle-filtering algorithms, the following steps are also necessary:
-   Define a corresponding [`State`](@ref) sub-type;
-   (optional) Define `Patter.map_init()` and `Patter.states_init()` methods for [`simulate_states_init()`](@ref) to simulate initial states;
-   Define a [`Patter.simulate_step()`](@ref) method (for [`Patter.simulate_move()`](@ref)) to update the state using a [`ModelMove`](@ref) instance (in [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref));
-   Define a [`Patter.logpdf_step()`](@ref) method (for [`Patter.logpdf_move()`](@ref)) to evaluate the probability density of movement from one state to another (in [`particle_smoother_two_filter()`](@ref));

"""
abstract type ModelMove end 

# Random walk in X and Y
struct ModelMoveXY{T, U, V, W} <: ModelMove
    # Environment
    map::T
    # Distribution of step lengths
    mobility::U
    dbn_length::V
    # Distribution of headings
    dbn_heading::W
end 

# Random walk in X, Y and Z
struct ModelMoveXYZ{T, U, V, W, X} <: ModelMove
    # Environment
    map::T
    # Step length
    mobility::U
    dbn_length::V
    # Heading
    dbn_heading::W
    # Depth
    dbn_z::X
end 

# Correlated random walk in XY
struct ModelMoveCXY{T, U, V, W} <: ModelMove
    # Environment
    map::T
    # Steplength 
    mobility::U
    dbn_length::V
    # Turning angle 
    # * This must be symmetrically centred around zero
    dbn_heading_delta::W
end 

# Correlated random walk in XYZ
struct ModelMoveCXYZ{T, U, V, W, X} <: ModelMove
    # Environment
    map::T
    # Step length
    mobility::U
    dbn_length::V
    # Turning angle
    # * This must be symmetrically centred around zero
    dbn_heading_delta::W
    # Change in depth 
    dbn_z_delta::X
end 


#########################
#########################
#### Simulate movement steps

"""
    simulate_step(state::State, model_move::ModelMove, t::Int64)

(Internal) Simulate a (tentative) step from one location ([`State`](@ref)) into a new location ([`State`](@ref)).

# Arguments

- `state`: A [`State`](@ref) instance that defines the animal's previous [`State`](@ref);
- `model_move`: A [`ModelMove`](@ref) instance;
- `t`: An integer that defines the time step;

# Details

[`Patter.simulate_step()`](@ref) is an internal generic function that simulates a new value for the animal's [`State`], that is, the animal's location (and other state components). Methods are provided for the built-in [`State`] and movement model ([`ModelMove`](@ref)) sub-types. For custom [`State`](@ref)s or [`ModelMove`](@ref) sub-types, corresponding methods must be provided. Internally, [`Patter.simulate_step()`](@ref) is wrapped by [`Patter.simulate_move()`](@ref), which implements[ `simulate_step()`](@ref) iteratively until a valid [`State`](@ref) is simulated (see [`is_valid()`](@ref)). 

# Returns

- A [`State`](@ref) instance;

# See also

* [`State`](@ref) and [`ModelMove`](@ref) for [`State`](@ref) and movement model sub-types;
* [`Patter.simulate_step()`](@ref) to simulate a new [`State`](@ref);
* [`Patter.is_valid()`](@ref) to determine whether or not a simulated state is valid;
* [`Patter.simulate_move()`](@ref) to simulate states iteratively until a valid state is found;
* [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref) for the front-end functions that use these routines to simulate animal movement paths;

"""
function simulate_step end

# RW in XY
function simulate_step(state::StateXY, model_move::ModelMoveXY, t::Int64)
    length    = rand(model_move.dbn_length)
    heading   = rand(model_move.dbn_heading)
    x         = state.x + length * cos(heading)
    y         = state.y + length * sin(heading)
    map_value = extract(model_move.map, x, y)
    StateXY(map_value, x, y)
end 

# RW in XYZ
function simulate_step(state::StateXYZ, model_move::ModelMoveXYZ, t::Int64)
    length    = rand(model_move.dbn_length)
    heading   = rand(model_move.dbn_heading)
    x         = state.x + length * cos(heading)
    y         = state.y + length * sin(heading)
    z         = rand(model_move.dbn_z)
    map_value = extract(model_move.map, x, y)
    StateXYZ(map_value, x, y, z)
end 

# CRW in XY
function simulate_step(state::StateCXY, model_move::ModelMoveCXY, t::Int64)
    length  = rand(model_move.dbn_length)
    heading = state.heading + rand(model_move.dbn_heading_delta)
    x       = state.x + length * cos(heading)
    y       = state.y + length * sin(heading)
    map_value = extract(model_move.map, x, y)
    StateCXY(map_value, x, y, heading)
end 

# CRW in XYZ 
function simulate_step(state::StateCXYZ, model_move::ModelMoveCXYZ, t::Int64)
    length  = rand(model_move.dbn_length)
    heading = state.heading + rand(model_move.dbn_heading_delta)
    x       = state.x + length * cos(heading)
    y       = state.y + length * sin(heading)
    z       = state.z + rand(model_move.dbn_z_delta)
    map_value = extract(model_move.map, x, y)
    StateCXYZ(map_value, x, y, z, heading)
end 


"""
    simulate_move(state::State, model_move::ModelMove, t::Int64, n_trial::Real)

(Internal) Simulate movement from one location ([`State`](@ref)) into a new location ([`State`](@ref)).

# Arguments

- `state`: A [`State`](@ref) instance that defines the animal's previous [`State`](@ref);
- `model_move`: A [`ModelMove`](@ref) instance;
- `t`: An integer that defines the time step;
- `n_trial`: A number that defines the number of attempts to simulate a valid [`State`](@ref);

# Details

[`Patter.simulate_move()`](@ref) is an internal function that uses a [`Patter.simulate_step()`](@ref) method to simulate new state proposals iteratively until a valid [`State`](@ref) is generated or `n_trial` is reached (see [`is_valid()`](@ref)). For custom [`State`](@ref)s or [`ModelMove`](@ref) sub-types, corresponding [`Patter.simulate_step()`](@ref) methods must be provided for this function. [`Patter.simulate_move()`](@ref) is used to simulate movement paths (e.g., in [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref)).

# Returns

- A `Tuple` that comprises the simulated [`State`](@ref) instance and the (log) weight; i.e., 
    * `([State], 0.0)` if [`State`](@ref) is valid;
    * `([State], -Inf)` if [`State`](@ref) is invalid;

# See also

* [`State`](@ref) and [`ModelMove`](@ref) for [`State`](@ref) and movement model sub-types;
* [`Patter.simulate_step()`](@ref) to simulate a new [`State`](@ref);
* [`Patter.is_valid()`](@ref) to determine whether or not a simulated state is valid;
* [`Patter.simulate_move()`](@ref) to simulate states iteratively until a valid state is found;
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
    logpdf_step(state_from::State, state_to::State, model_move::ModelMove, length, heading)
    
(Internal) Evaluate the (unnormalised) log probability of an (unrestricted) movement step. 

# Arguments

- `state_from`: A [`State`](@ref) instance that defines a [`State`](@ref) from which the animal moved;
- `state_to`: A [`State`](@ref) instance that defines a [`State`](@ref) into which the animal moved;
- `model_move`: A [`ModelMove`](@ref) instance;
- `t`: An integer that defines the time step;
- `length`: A float that defines the step length (i.e., the Euclidean distance between `state_from` (`x`, `y`) and `state_to` (`x`, `y`));
- `heading`: A float that defines the angle (in polar coordinates) between `state_from` (`x`, `y`) and `state_to` (`x`, `y`);

# Details

[`Patter.logpdf_step()`](@ref) is an internal generic function that evaluates the (unnormalised) log probability of an (unrestricted) movement step between two [`State`](@ref)(s) (i.e., locations). Methods are provided for the built-in [`State`](@ref) and [`ModelMove`](@ref) sub-types, but need to be provided for custom sub-types. Internally, [`Patter.logpdf_step()`](@ref) is wrapped by [`Patter.logpdf_move()`](@ref), which evaluates the log probability of movement between two [`State`](@ref)s, accounting for restrictions to movement; that is, [`Patter.logpdf_move()`](@ref) evaluates `logpdf_step(state_from, state_to, model_move, length, heading) + log(abs(determinate)) - log(Z)` where `Z` is the normalisation constant. This is required for particle smoothing (see [`particle_smoother_two_filter()`](@ref)).

# Returns

- A number (log probability); 

# See also 

* [`State`](@ref) and [`ModelMove`](@ref) for [`State`](@ref) and movement model sub-types;
* [`Patter.simulate_step()`](@ref) and [`Patter.simulate_move()`](@ref) to simulate new [`State`](@ref)s;
* [`Patter.logpdf_step()`](@ref) and [`Patter.logpdf_move()`](@ref) to evaluate the log-probability of movement between two locations;
* [`Patter.logpdf_move_normalisation()`](@ref) for estimation of the normalisation constant;
* [`particle_smoother_two_filter()`](@ref) for the front-end function that uses these routines for particle smoothing;

"""
function logpdf_step end 

function logpdf_step(state_from::StateXY, state_to::StateXY, model_move::ModelMoveXY, t::Int64, length::Float64, heading::Float64) 
    logpdf(model_move.dbn_length, length) + logpdf(model_move.dbn_heading, heading)
end 

function logpdf_step(state_from::StateXYZ, state_to::StateXYZ, model_move::ModelMoveXYZ, t::Int64, length::Float64, heading::Float64) 
    logpdf(model_move.dbn_length, length) + logpdf(model_move.dbn_heading, heading) + logpdf(model_move.dbn_z, state_to.z)
end 

function logpdf_step(state_from::StateCXY, state_to::StateCXY, model_move::ModelMoveCXY, t::Int64, length::Float64, heading::Float64) 
    # Compute turning angle 
    heading_delta = abs_angle_difference(heading, state_from.heading)
    # Sum up logpdfs
    logpdf(model_move.dbn_length, length) + logpdf(model_move.dbn_heading_delta, heading_delta)
end 

function logpdf_step(state_from::StateCXYZ, state_to::StateCXYZ, model_move::ModelMoveCXYZ, t::Int64, length::Float64, heading::Float64) 
    # Compute change in depth 
    z_delta = state_to.z - state_from.z
    # Compute turning angle 
    heading_delta = abs_angle_difference(heading, state_from.heading)
    # Sum up logpdfs
    logpdf(model_move.dbn_length, length) + logpdf(model_move.dbn_heading_delta, heading_delta) + logpdf(model_move.dbn_z_delta, z_delta)
end 


"""
    logpdf_move(state_from::State, state_to::State, state_zdim::Bool, 
                model_move::ModelMove, 
                t::Int, 
                vmap::Union{GeoArray, Nothing}, 
                n_sim::Int,
                cache::Union{Dict, Nothing})
    
(Internal) Evaluate the log probability of a movement between two [`State`](@ref)s (`state_from` and `state_to`). 

# Arguments

- `state_from`: A [`State`](@ref) instance that defines a [`State`](@ref) from which the animal moved;
- `state_to`: A [`State`](@ref) instance that defines a [`State`](@ref) into which the animal moved;
- `state_zdim`: A `Boolian` that defines whether or not `state_from` and `state_to` contain a `z` (depth) dimension;
- `model_move`: A [`ModelMove`](@ref) instance;
- `t`: An integer that defines the time step;
- `vmap`: (optional) A `GeoArray` that maps the region within which movements from `state_from` are always legal. Valid regions must equal 1. `vmap` can be provided for 'horizontal' movement models (e.g., if `state_from` and `state_to` are `StateXY` instances (see [`State`](@ref)));
- `n_sim`: An integer that defines the number of Monte Carlo simulations (used to approximate the normalisation constant);
- `cache`: (optional) A Dict of normalisation constants including `state_from`;

# Details

[`Patter.logpdf_move()`](@ref) is an internal function that evaluates the log probability of a movement step between two [`State`](@ref)(s) (i.e., locations). This function wraps [`Patter.logpdf_step()`](@ref), accounting for accounting for restrictions to movement; that is, [`Patter.logpdf_move()`](@ref) evaluates `logpdf_step(state_from, state_to, model_move, t, length, heading) + log(abs(determinate)) - log(Z)` where `Z` is the normalisation constant. If `model_move` is 'horizontal (e.g., `state_from` and `state_to` are two-dimensional, `StateXY` instances), a 'validity map' (`vmap`) can be provided. This is a `GeoArray` that define the regions within which movements between two locations are always legal. In the case of an aquatic animal, this is the region of the study area that is the sea, shrunk by `state_from.mobility`. In this instance, the normalisation constant is simply `log(1.0)`. Otherwise, a Monte Carlo simulation of `n_sim` iterations is required to approximate the normalisation constant, accounting for invalid movements, which is more expensive (see [`logpdf_move_normalisation()`](@ref)). Set `n_sim = 0` and `cache = nothing` for unrestricted models (i.e., if [`particle_filter()`](@ref) was implemented with `n_move = 1`). [`Patter.logpdf_move()`](@ref) is used for particle smoothing (see [`particle_smoother_two_filter()`](@ref)).

# Returns

- A number (log probability); 

# See also 

* [`State`](@ref) and [`ModelMove`](@ref) for [`State`](@ref) and movement model sub-types;
* [`Patter.simulate_step()`](@ref) and [`Patter.simulate_move()`](@ref) to simulate new [`State`](@ref)s;
* [`Patter.logpdf_step()`](@ref) and [`Patter.logpdf_move()`](@ref) to evaluate the log-probability of movement between two locations;
* [`Patter.logpdf_move_normalisation()`](@ref) for estimation of the normalisation constant;
* [`particle_smoother_two_filter()`](@ref) for the front-end function that uses these routines for particle smoothing;

"""
function logpdf_move(state_from::State, state_to::State, state_zdim::Bool, 
                     model_move::ModelMove, t::Int, vmap::Union{GeoArray, Nothing}, n_sim::Int, 
                     cache::Union{Dict{<:State, Float64}, Nothing}) 

    #### Validate state 
    # (This check is cheap)
    if !state_is_valid(state_to, state_zdim)
         return -Inf
    end

    #### Evaluate step lengths and headings 
    x = state_to.x - state_from.x
    y = state_to.y - state_from.y
    length, heading = cartesian_to_polar(x, y)
    # When locations exceed model_move.mobility, we set -Inf density for speed
    # * Consider NearestNeighbors.jl to build a kd tree:
    # * - Provide locations
    # * - Get index of particles within mobility
    # * - This may improve speed by eliminating the need to iterative over all particles in the smoother 
    if length > model_move.mobility
        return -Inf
    end
    # Calculate log(abs(determinate))
    log_det = -0.5 * log(x^2 + y^2)
    
    #### Compute normalisation constant
    if !isnothing(cache)
        # Lookup normalisation constant from cache
        log_z = getindex(cache, state_from)
    else 
        # Compute normalisation constant on the fly
        # * E.g., For movement models that depend on both `state_from` and `t`
        log_z = logpdf_move_normalisation(state_from, state_zdim, model_move, t, vmap, n_sim)
    end 

    #### Evaluate density 
    logpdf_step(state_from, state_to, model_move, t, length, heading) + log_det - log_z

end 


"""
    logpdf_move_normalisation(state::State, state_zdim::Bool, 
                              model_move::ModelMove, t::Int, 
                              vmap::Union{GeoArray, Nothing},
                              n_sim::Int)

(Internal) Approximate the (log) normalisation constant for the (log) probability density of movement from one [`State`](@ref) (location) into another. 

# Arguments

- `state`: A [`State`](@ref) instance that defines a [`State`](@ref) from which the animal moved;
- `state_zdim`: A `Boolian` that defines whether or not `state` contains a `z` (depth) dimension;
- `model_move`: A [`ModelMove`](@ref) instance;
- `t`: An integer that defines the time step;
- `vmap`: (optional) A `GeoArray` that maps the region within which movements from `state` are always legal. Valid regions must equal 1. `vmap` can be provided for 'horizontal' movement models (e.g., if `state` is a `StateXY`);
- `n_sim`: An integer that defines the number of Monte Carlo simulations. Set `n_sim = 0` if `n_move = 1` in [`particle_filter()`](@ref);

# Details

This internal function computes the normalisation constant for the (log) probability of movement from one [`State`](@ref) (`state`) into another. If `n_move = 1` in [`particle_filter()`](@ref), set `n_sim = 0` to return `log(1.0)`. Otherwise, we need to account for the truncation of the movement model by land. If `model_move` is 'horizontal' (e.g., `state` is a two-dimensional, `StateXY` instance), a 'validity map' (`vmap`) can be provided. This is a `GeoArray` that define the regions within which movements from that `state` are always legal. In the case of an aquatic animal, this is the region of the study area that is the sea, shrunk by `state.mobility`. In this instance, the normalisation constant is simply `log(1.0)`. Otherwise, a Monte Carlo simulation of `n_sim` iterations is used to estimate the normalisation constant. A Beta(1, 1) prior is used to correct for simulations that fail to generate valid move from `state`. This function is used by [`Patter.logpdf_move()`](@ref) to evaluate the (log) probability of movement between two states, which is required for particle smoothing (see [`particle_smoother_two_filter()`](@ref)).

# Returns 

- A number (the log normalisation constant); 

# See also

* [`State`](@ref) and [`ModelMove`](@ref) for [`State`](@ref) and movement model sub-types;
* [`Patter.simulate_step()`](@ref) and [`Patter.simulate_move()`](@ref) to simulate new [`State`](@ref)s;
* [`Patter.logpdf_step()`](@ref) and [`Patter.logpdf_move()`](@ref) to evaluate the log-probability of movement between two locations;
* [`Patter.logpdf_move_normalisation()`](@ref) for estimation of the normalisation constant;
* [`particle_smoother_two_filter()`](@ref) for the front-end function that uses these routines for particle smoothing;

"""
function logpdf_move_normalisation(state::State, state_zdim::Bool, model_move::ModelMove, t::Int, vmap::Union{GeoArray, Nothing}, n_sim::Int)

    # (A) Set normalisation constant to zero for n_sim = 0
    # * This corresponds to n_move = 1 in particle_filter()
    if n_sim == 0
        return 0.0
    end

    # (B) Extract normalisation constant from vmap, if possible 
    if !isnothing(vmap) && isone(extract(vmap, state.x, state.y))
        # (A) Use vmap 
        # * vmap may be supplied for 'horizontal' (2D) movement models 
        # * vmap elements = 0 or 1 
        # * vmap = 1.0 indicates that all moves from that point are valid (land in water)
        log_z = 0.0 

    else 
        # (C) Compute normalisation constant via simulation
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
        log_z = log((k + 1) / (n_sim + 2))
    end 

    return log_z

end

# Precompute normalisation constants for all unique states 
function logpdf_move_normalisations(states::Matrix{<:State}, model_move::ModelMove, vmap::Union{GeoArray, Nothing}, n_sim::Int)

    # Identify dimension of input state
    state_zdim = hasfield(typeof(states[1]), :z)
 
    # Compute normalisation constants
    # * If the cache is used, `t` is assumed irrelevant and simply set to 1 here
    # * This operation is faster if multi-threaded
    unique_states = unique(states)
    constants = Vector{Pair{eltype(unique_states), Float64}}(undef, length(unique_states))
    @threads for i in eachindex(unique_states)
        constants[i] = unique_states[i] => logpdf_move_normalisation(unique_states[i], state_zdim, model_move, 1, vmap, n_sim)
    end
    
    # Move constants into dictionary
    cache = Dict(constants)

    # Return cache 
    return cache 

end 