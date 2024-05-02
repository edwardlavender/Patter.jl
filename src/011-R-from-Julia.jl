using DataFrames

"""
# R from Julia

A collection of functions that facilitate the translation of inputs from `Julia` into `R`. 

# Details
* [`r_get_states`] translates a State matrix into a `DataFrame` that can be passed to `R`. In the input matrix, each row is a particle and each column is a time step. 

# Returns 
A long-format `DataFrame`, with columns for `path_id`, `timestep` and each state dimension.
"""
function r_get_states(state::Matrix)
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
    # Coerce to dataframe
    fields = (:path_id, :timestep, fields...)
    DataFrame(values, collect(fields))
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