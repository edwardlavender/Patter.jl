using Dates

export assemble_yobs

#########################
#########################
#### Dictionary helpers

# Initialise the entry of a dictionary given observations 
# * This method is used for data preparation
# * The user provides a vector of dictionaries that contain observations and a corresponding character vector of models
# * We extract the observation type & model types for each vector element & formulate the dictionary entry
function dict_initialise_entry(datasets::Vector, models::Vector{String})
    obs_type = [typeof(dataset.obs[1]) for dataset in datasets]
    mod_type = [eval(Symbol(model)) for model in models]
    tuples   = [Tuple{obs_type[i], mod_type[i]} for i in eachindex(obs_type)]
    Vector{Union{tuples...}}
end 

# Initialise a dictionary with entries of a given type 
# * `entry` is specified via dict_initialise_entry()
function dict_initalise(entry) 
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
-   models: A Vector of Strings that define the ModelObs subtype for each dataset. 

# Details

The function iterates over datasets and models and creates a typed dictionary of timestamps. Each time step contains a vector of Tuples, with one element for each sensor that recorded an observation at that time stamp. Each element is a tuple that defines the observation and the model parameters (that is, a `ModelObs` instance).

"""
function assemble_yobs(datasets, models) 

    # Initialise empty dictionary 
    # ydict = Dict{DateTime, Vector{Union{Tuple{Float64, ModelObsDepthUniform}, Tuple{Int64, ModelObsAcousticLogisTrunc}}}}()
    entry = dict_initialise_entry(datasets, models)
    ydict = dict_initalise(entry)
    entry = entry[]

    # Iterate over datasets/models
    for (dataset, model) in zip(datasets, models)
    
        # Define dataset parameters and ModelObs subtype
        parameters      = select(dataset, Not([:timestamp, :sensor_id, :obs]))
        # strings         = any(col -> isa(first(col), String), eachcol(dataset))
        ModelObsSensor  = eval(Symbol(model))
    
        # Iterate over time steps & add observations & ModelObs objects
        for (timestamp, obs, row) in zip(dataset.timestamp, dataset.obs, Tables.namedtupleiterator(parameters))
            if !haskey(ydict, timestamp)
                # ydict[timestamp] = Vector{Union{Tuple{Float64, ModelObsDepthUniform}, Tuple{Int64, ModelObsAcousticLogisTrunc}}}() 
                ydict[timestamp] = entry
            end 
            # Build ModelObs object
            # values = row
            # if (strings)
            #    row = map(eval_string_as_symbol, row)
            # end            
            # Update dictionary 
            push!(ydict[timestamp], (obs, ModelObsSensor(row...)))
        end 
    
    end 
    
    ydict

end