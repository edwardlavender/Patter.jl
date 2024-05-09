using DataFrames
using OrderedCollections

"""
# R from Julia

A collection of functions that facilitate the translation of inputs from `Julia` into `R`. 

# Details
* [`r_get_states`] translates a State matrix into a `DataFrame` that can be passed to `R`. In the input matrix, each row is a particle and each column is a time step. 
* [`r_get_dataset`] translates a Dictionary of observations into a Vector of DataFrames that can be passed to `R`.

# Returns 
A long-format `DataFrame`, with columns for `path_id`, `timestep` and each state dimension.
"""
function r_get_states(state::Matrix, timesteps::Vector = collect(1:size(state, 2)))
    # Initialise empty matrix
    fields = fieldnames(typeof(state[1]))
    values = Matrix{Float64}(undef, prod(size(state)), length(fields) + 2)
    # Define path ID & time step columns
    np = size(state, 1)
    nt = size(state, 2)
    values[:, 1] = repeat(1:np, inner = nt)
    values[:, 2] = repeat(1:nt, outer = np)
    # Populate matrix
    for i in 1:size(values, 1)
        for j in eachindex(fields)
            values[i, j + 2] = getfield(state[Int(values[i, 1]), Int(values[i, 2])], fields[j])
        end 
    end 
    # Replace column indices with timesteps 
    values[:, 2] = repeat(minimum(timesteps):maximum(timesteps), outer = np)
    # Coerce to dataframe
    fields = (:path_id, :timestep, fields...)
    df = DataFrame(values, collect(fields))
    df.path_id = Int.(df.path_id)
    df.timestep = Int.(df.timestep)
    df
end

# Examples:

# Define state matrix: 
# * Two rows: two particles
# * Three columns: three time steps 
# state = [StateXY(0.0, 1.0, 2.0)  StateXY(0.0, 3.0, 4.0)  StateXY(0.0, 5.0, 6.0);
# StateXY(0.0, 7.0, 8.0) StateXY(0.0, 9.0, 10.0)    StateXY(0.0, 11.0, 2.0)]
# r_get_states(state)

# Create a big matrix of StateXY objects
# np = 1000
# nt = 20000
# state = [StateXY(rand(), rand(), rand()) for _ in 1:np, _ in 1:nt]
# r_get_states(state)

#### Convert struct to dict
function struct_to_dict(s)
    return OrderedCollections.OrderedDict(key => getfield(s, key) for key in propertynames(s))
end

#### Formulate dictionaries to hold the observation for a selected time stamp
function dict_obs(timestamp, obs, sensor)
    OrderedCollections.OrderedDict(:timestamp => timestamp, :obs => obs, struct_to_dict(sensor)...)
end 

#### Extract dataset(s) from yobs for a selected model type
function r_get_dataset(yobs::Dict, model_type::Type{<: ModelObs})

    # Initialise a Vector of DataFrames 
    # * One dataframe for each dictionary element (i.e., simulated path)
    # * Note that output = [] is required for correct translation via JuliaCall 
    output = [] # Vector{DataFrame}()

    # Build the vector of dataframes
    for key in sort(collect(keys(yobs)))
        # Initialise a Vector of dictionaries
        # * Each dictionary will contain, for a specified time stamp, the observation & essential model parameters 
        dicts = Vector{OrderedDict}()
        for (timestamp, observations) in yobs[key]
            for (obs, sensor) in observations
                if sensor isa model_type
                    push!(dicts, dict_obs(timestamp, obs, sensor))
                end
            end 
        end
        # Convert the dictionaries to a DataFrame
        df = DataFrame(dicts)
        # df.key .= key 
        push!(output, df)
    end 

    output

end 

function r_get_particles(particles::NamedTuple)
    # Collate information 
    states      = r_get_states(particles.state, particles.timesteps)
    diagnostics = DataFrame(timestep = particles.timesteps, 
                            timestamp = particles.timestamps, 
                            ess = particles.ess, 
                            maxlp = particles.maxlp)
    # Return outputs 
    (
        states     = states, 
        diagnostics = diagnostics, 
        convergence = particles.convergence
    )
end 