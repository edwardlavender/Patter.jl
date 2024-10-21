using Dates

export simulate_path_walk, simulate_yobs

"""
    simulate_path_walk(; xinit = Vector, model_move::ModelMove, timeline::Vector{DateTime})

Simulate discrete-time movement path(s) from a `Vector` of initial [`State`](@ref)s and a random-walk movement model. 

# Arguments (keywords)

- `xinit`: A Vector of initial [`State`](@ref)s instances;
- `model_move`: A [`ModelMove`](@ref) instance;
- `timeline`: A `Vector{DateTime}` of ordered, regularly spaced time stamps that defines the time steps for the simulation;

# Details

[`State`](@ref) refers to the (`x`, `y`) location of an animal (alongside additional state components, if applicable). To simulate initial states, use [`simulate_states_init()`](@ref). For each initial state, [`simulate_path_walk()`](@ref) simulates a sequence of [`State`](@ref)s (i.e., a movement path) of `length(timeline)` steps using the movement model (`model_move`). The simulation of movement from one [`State`](@ref) into another is implemented by the internal function [`Patter.simulate_move()`](@ref), which in turn wraps [`Patter.simulate_step()`](@ref). At each time step, [`Patter.simulate_move()`](@ref) implements [`Patter.simulate_step()`](@ref) iteratively until a valid movement is identified (see [`is_valid()`](@ref)). [`Patter.simulate_step()`](@ref) is a generic function. Methods are implemented for the built-in [`State`](@ref) and [`ModelMove`](@ref) sub-types but custom sub-types require a corresponding [`Patter.simulate_step()`](@ref) method. 

# Returns
- A `matrix` of [`State`](@ref)s:
    - Each row represents a simulated path;
    - Each column represents a time step along `timeline`;

# See also

"""
function simulate_path_walk(; xinit = Vector, model_move::ModelMove, timeline::Vector{DateTime})

    # Define output objects
    np = length(xinit)
    nt = length(timeline)
    xout = Matrix{eltype(xinit)}(undef, np, nt);
    xout[:, 1] .= xinit
    
    # Run simulation
    for t in 2:nt
        for i in 1:np
            xi, lwi = simulate_move(xout[i, t - 1], model_move, t, 1_000_000)
            if isinf(lwi)
                error("`simulate_path_walk()` failed to simulate a valid step at time $t after 1 million trials. It is likely that `model_move` (which includes the study domain [map] and the movement model) has been incorrectly specified.")
            end 
            xout[i, t] = xi
        end
    end 

    xout

end


"""
    simulate_yobs(; paths::Matrix, model_obs::Vector{ModelObs}, timeline::Vector{DateTime})

For a series of simulated paths, simulate a dictionary of observations. 

# Arguments

- `paths`: A `Matrix` of simulated paths from [`simulate_path_walk()`](@ref);
- `model_obs`: A Vector of [`ModelObs`](@ref) instances;
- `timeline`: A `Vector{DateTime}` of ordered, regularly spaced time stamps that defines the time steps for the simulation;

# Details

The function expects a `Matrix` of simulated paths (see [`simulate_path_walk()`](@ref)). For each simulated path, the function iterates over each step in `timeline` and simulates observations using the `Vector` of observation models. Observations are simulated by the internal generic [`simulate_obs()`](@ref) via `simulate_obs(State, model, t)`, where `t` is the time step. Methods are provided for the built-in [`State`](@ref)s and [`ModelObs`](@ref) sub-types. For custom sub-types, a corresponding [`simulate_obs()`](@ref) method is required. Simulated observations can be used in the particle filter to reconstruct the underlying movements (see [`particle_filter()`](@ref)).

# Returns
- A `Dict`, with one entry for each path:
    - Each entry is a `Dict`, with one entry for each time stamp;
        - Each time stamp entry is a `Vector` of `Tuple`s, each comprising the simulated observation and the associated [`ModelObs`](@ref) instance (see also [`assemble_yobs()`](@ref));

# See also 

* [`State`](@ref) and [`ModelMove`](@ref) for [`State`](@ref) and movement model sub-types;
* [`Patter.simulate_step()`](@ref) and [`Patter.simulate_move()`](@ref) to simulate new [`State`](@ref)s;
* [`simulate_path_walk()`](@ref) to simulate animal movement paths (via [`ModelMove`](@ref));
* [`simulate_yobs()`](@ref) to simulate observations arising from simulated movements (via [`ModelObs`](@ref));

"""
function simulate_yobs(; paths::Matrix, model_obs::Vector{ModelObs}, timeline::Vector{DateTime})
    
    #### Initialise a set of dictionaries (one per path)
    # Initialise a set of dictionaries
    entry          = dict_initialise_entry(paths[1, 1], model_obs)
    yobs_by_path   = Dict{Int, Dict{DateTime, entry}}()
    # Define the initialisation object for each path entry
    path_dict_init = dict_initialise(entry)

    #### Iterate over paths & populate dictionaries 
    for i in 1:size(paths, 1)
        
        # Initialise dictionary for path
        path = paths[i, :]
        yobs_by_path[i] = deepcopy(path_dict_init)

        # Simulate observations 
        for (t, state) in enumerate(path)
            yobs_by_path[i][timeline[t]] = entry[]
            for model in model_obs
                push!(yobs_by_path[i][timeline[t]], (simulate_obs(state, model, t), model))
            end
        end

    end 

    yobs_by_path

end