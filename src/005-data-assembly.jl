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
    # Assemble `yobs` for the particle filter

Assemble a dictionary of observations (and associated model parameters) for the particle filter. 

# Arguments
-   datasets: A Vector of DataFrames, one for each data type. Each DataFrame must contain the following columns:
    - `timestamp`: A DateTime Vector of time stamps;
    - `sensor_id`: A vector of sensor IDs;
    - `obs`: The observation;
    - Additional columns required to construct ModelObs instances (that is, model parameters);
-   model_types: A Vector of ModelObs subtypes for each dataset. 

# Details

The function iterates over datasets and models and creates a typed dictionary of timestamps. Each time step contains a vector of Tuples, with one element for each sensor that recorded an observation at that time stamp. Each element is a tuple that defines the observation and the model parameters (that is, a `ModelObs` instance).

"""
function assemble_yobs(datasets::Vector, model_types::Vector{DataType}) 

    # Initialise empty dictionary 
    # yobs = Dict{DateTime, Vector{Union{Tuple{Float64, ModelObsDepthUniform}, Tuple{Int64, ModelObsAcousticLogisTrunc}}}}()
    entry = dict_initialise_entry(datasets, model_types)
    yobs = dict_initialise(entry)
    entry = entry[]

    # Iterate over datasets/models
    for (dataset, model_type) in zip(datasets, model_types)
    
        # Define dataset parameters and ModelObs subtype
        parameters      = select(dataset, Not([:timestamp, :obs]))
        # strings         = any(col -> isa(first(col), String), eachcol(dataset))
    
        # Iterate over time steps & add observations & ModelObs objects
        for (timestamp, obs, row) in zip(dataset.timestamp, dataset.obs, Tables.namedtupleiterator(parameters))
            if !haskey(yobs, timestamp)
                # yobs[timestamp] = Vector{Union{Tuple{Float64, ModelObsDepthUniform}, Tuple{Int64, ModelObsAcousticLogisTrunc}}}() 
                yobs[timestamp] = entry
            end 
            # Build ModelObs object
            # values = row
            # if (strings)
            #    row = map(eval_string_as_symbol, row)
            # end            
            # Update dictionary 
            push!(yobs[timestamp], (obs, model_type(row...)))
        end 
    
    end 
    
    yobs

end