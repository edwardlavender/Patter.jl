using Dates

export simulate_path_walk, simulate_yobs

"""
# Simulate movement paths

This function simulates discrete-time movement path(s) from a vector of initial states and random-walk movement model. 

# Arguments
- `xinit`: A Vector of [`State`](@ref)s that defines the initial state(s) for the simulation;
- `move`: A [`ModelMove`](@ref) instance;
- `nt`: An integer that defines the number of time steps;

# Details

For each initial state, a movement path is simulated. 

# Returns
A matrix:
- Rows: paths
- Columns: time steps
"""
function simulate_path_walk(; xinit = Vector, move::ModelMove, timeline::Vector{DateTime})

    # Define output objects
    np = length(xinit)
    nt = length(timeline)
    xout = Matrix{eltype(xinit)}(undef, np, nt);
    xout[:, 1] .= xinit
    
    # Run simulation
    for t in 2:nt
        for i in 1:np
            xout[i, t] = simulate_move(xout[i, t - 1], move, t, Inf)[1]
        end
    end 

    xout

end


"""
# Simulate observations 

For each simulated path, simulate a dictionary of observations. 

# Arguments
- `paths`: A Matrix of simulated paths; from [`simulate_path_walk()`](@ref);
- `models`: A Vector of `ModelObs` instances;
- `timeline`: A Vector{DateTime} of time stamps;

# Returns
- A dictionary, with one entry for each path;
- Each entry is a dictionary, with one entry for each time stamp;
- Each time stamp entry is a Vector of Tuples, each comprising the observation and the associated `ModelObs` instance;

"""
function simulate_yobs(; paths::Matrix, models::Vector{ModelObs}, timeline::Vector{DateTime})
    
    #### Initialise a set of dictionaries (one per path)
    # Initialise a set of dictionaries
    entry          = dict_initialise_entry(paths[1, 1], models)
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
            for model in models
                push!(yobs_by_path[i][timeline[t]], (simulate_obs(state, model, t), model))
            end
        end

    end 

    yobs_by_path

end