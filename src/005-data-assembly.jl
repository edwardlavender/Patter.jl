using Dates

export assemble_yobs

#########################
#########################
#### Dictionary helpers

# Initialise the entry of a dictionary given a State
# * This method is used in simulations
# * The user specifies the State instance and the model instances
# * We use the State and the models to determine the observation type & formulate the dictionary entry
function dict_initialise_entry(state::State, models::Vector{ModelObs})
    obs_types   = [typeof(simulate_obs(state, model, 1)) for model in models]
    model_types = [typeof(model) for model in models]
    tuples      = [Tuple{obs_types[i], model_types[i]} for i in eachindex(obs_types)]
    Vector{Union{tuples...}}
end

# Initialise the entry of a dictionary given observations 
# * This method is used for data preparation
# * The user provides a vector of dictionaries that contain observations and a corresponding vector of model types
# * We extract the observation type & model types for each vector element & formulate the dictionary entry
function dict_initialise_entry(datasets::Vector, model_types::Vector{DataType})
    obs_types = [typeof(dataset.obs[1]) for dataset in datasets]
    tuples   = [Tuple{obs_types[i], model_types[i]} for i in eachindex(obs_types)]
    Vector{Union{tuples...}}
end 

# Initialise a dictionary with entries of a given type 
# * `entry` is specified via dict_initialise_entry()
function dict_initialise(entry) 
    Dict{DateTime, entry}()
end 


#########################
#########################
#### Data assembly

"""
    assemble_yobs(; datasets::Vector, model_obs_types::Vector{DataType})

Assemble a dictionary of observations (and associated model parameters) for the particle filter ([`particle_filter()`](@ref)). 

# Arguments (keywords)

-   `datasets`: A `Vector` of `DataFrame`s, one for each data type. Each `DataFrame` must contain the following columns:
    - `timestamp`: A `DateTime` `Vector` of time stamps;
    - `sensor_id`: A `Vector` of sensor IDs;
    - `obs`: A `Vector` of observations;
    - Additional columns required to construct [`ModelObs`](@ref) instances (that is, model parameters);
-   `model_obs_types`: A `Vector` of [`ModelObs`](@ref) sub-types (one for each `dataset`);

# Details

The function iterates over animal-tracking datasets (for example, acoustic and archival [depth] time series _for a particular individual_) and corresponding observation model types and creates a typed dictionary of time stamps for the particle filter (see [`particle_filter()`](@ref)'s `yobs` argument). Each time step contains a `Vector` of `Tuple`s, with one element for each sensor that recorded an observation at that time stamp. For example, each element might correspond to an acoustic receiver and/or the individual's archival tag. Each element is a `Tuple` that defines the observation and the corresponding observation model parameters (that is, a [`ModelObs`](@ref) instance). 

# Returns

- A `Dict`:
    - Each element corresponds to a `timestamp`:
        - Each timestamped element contains a `Vector` of `Tuples` (one for each observation):
            - Each `Tuple` contains an observation and the corresponding [`ModelObs`](@ref) instance;

# See also 

* [`assemble_yobs()`](@ref) to assemble real-world datasets for the particle filter;
* [`simulate_yobs()`](@ref) to simulate observations for the particle filter;
* [`particle_filter()`](@ref) to implement the particle filter;

"""
function assemble_yobs(; datasets::Vector, model_obs_types::Vector{DataType}) 

    # TO DO
    # * Review the speed of this function

    # Initialise empty dictionary 
    # yobs = Dict{DateTime, Vector{Union{Tuple{Float64, ModelObsDepthUniform}, Tuple{Int64, ModelObsAcousticLogisTrunc}}}}()
    entry = dict_initialise_entry(datasets, model_obs_types)
    yobs = dict_initialise(entry)
    entry = entry[]

    # Iterate over datasets/models
    for (dataset, model_obs_type) in zip(datasets, model_obs_types)
    
        # Define dataset parameters
        parameters = select(dataset, Not([:timestamp, :obs]))
        # Enforce a column order that matches the components of model_obs_type
        parameters = parameters[:, [col for col in fieldnames(model_obs_type)]] 
    
        # Iterate over time steps & add observations & ModelObs objects
        for (timestamp, obs, row) in zip(dataset.timestamp, dataset.obs, Tables.namedtupleiterator(parameters))
            if !haskey(yobs, timestamp)
                # yobs[timestamp] = Vector{Union{Tuple{Float64, ModelObsDepthUniform}, Tuple{Int64, ModelObsAcousticLogisTrunc}}}() 
                yobs[timestamp] = entry
            end  
            # Update dictionary 
            push!(yobs[timestamp], (obs, model_obs_type(row...)))
        end 
    
    end 
    
    yobs

end