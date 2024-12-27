using DataFrames
using OrderedCollections

# """
# # `R` from `Julia`

# A collection of internal functions that facilitate the translation of `Julia` objects into `R`. 

# # Details

# * [`Patter.r_get_dataset()`](@ref) translates a Dictionary of observations into a `Vector` of `DataFrame`s that can be passed to `R`.
# * `Patter.r_get_states()` translates a `Matrix` of [`State`](@ref)s into a `DataFrame` that can be passed to `R`. In the input `Matrix`, each row is a particle and each column is a time step. 
# * [`Patter.r_get_particles()`](@ref) wraps `Patter.r_get_states()` and translates particle outputs (from [`particle_filter()`](@ref) and [`two_filter_smoother()`](@ref)) into a `NamedTuple` for `R`.

# These functions are [`State`](@ref) and model agnostic; that is, they work irrespective of the input [`State`](@ref) and model sub-types. Custom methods are not required to handle novel sub-types. 

# # Returns 

# * [`Patter.r_get_dataset()`](@ref) returns a `Vector` of `DataFrame`s, with columns for `timestamp`, `obs` and the observation model parameters;
# * `Patter.r_get_states()` returns a long-format `DataFrame`, with columns for `path_id`, `timestep` and each state dimension;
# * [`Patter.r_get_particles()`](@ref) returns a `NamedTuple` of particle information, including:
#     - `states`: A `DataFrame` of [`State`](@ref) dimensions (from `Patter.r_get_states()`);
#     - `diagnostics`: A `DataFrame` of algorithm diagnostics, including `timestep`, `timestamp`, `ess` and `maxlp` columns;
#     - `convergence`: A `Boolian` that defines algorithm convergence;

# """
# function r_get end 


# Convert a struct to an Ordered Dict
function struct_to_dict(s)
    return OrderedCollections.OrderedDict(key => getfield(s, key) for key in propertynames(s))
end


# Formulate dictionaries to hold the observation for a selected time stamp
function dict_obs(timestamp, obs, sensor)
    OrderedCollections.OrderedDict(:timestamp => timestamp, :obs => obs, struct_to_dict(sensor)...)
end 


# Extract dataset(s) from yobs for a selected model type
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
        sort!(df, [:timestamp, :sensor_id])
        # df.key .= key 
        push!(output, df)
    end 

    output

end 


# Get a DataFrame of States
function r_get_states(states::Matrix{<:State}, 
                      timesteps::Vector{Int},
                      timestamps::Vector{Dates.DateTime})
    # Initialise empty matrix
    if (length(timesteps) != length(timestamps))
        error("`length(timesteps)` and `length(timestamps)` should be identical.")
    end 
    fields = fieldnames(typeof(states[1]))
    values = Matrix{Float64}(undef, prod(size(states)), length(fields) + 2)
    # Define path ID & time step columns
    np = size(states, 1)
    nt = size(states, 2)
    values[:, 1] = repeat(1:np, inner = nt)
    values[:, 2] = repeat(1:nt, outer = np)
    # Populate matrix
    for i in 1:size(values, 1)
        for j in eachindex(fields)
            values[i, j + 2] = getfield(states[Int(values[i, 1]), Int(values[i, 2])], fields[j])
        end 
    end 
    # Replace column indices with timesteps 
    values[:, 2] = repeat(minimum(timesteps):maximum(timesteps), outer = np)
    # Coerce to dataframe
    df = DataFrame(values, collect((:path_id, :timestep, fields...)))
    df.path_id = Int.(df.path_id)
    df.timestep = Int.(df.timestep)
    df.timestamp = timestamps[indexin(df.timestep, timesteps)]
    select!(df, collect((:path_id, :timestep, :timestamp, fields...)))
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


# Get a Tuple of particle information (states, diagnostics, convergence)
function r_get_particles(particles::Particles)
    (
        states      = r_get_states(particles.states, 
                                   particles.diagnostics.timestep, 
                                   particles.diagnostics.timestamp), 
        diagnostics = particles.diagnostics, 
        callstats   = particles.callstats
    )
end 