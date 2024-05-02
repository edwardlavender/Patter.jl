using DataFrames

# Get state values from a State matrix > DataFrame
# (rows: particles, columns: time)
function rget_state_df(state::Matrix)
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
# state = [StateXY(1.0, 2.0)  StateXY(3.0, 4.0)  StateXY(5.0, 6.0);
# StateXY(7.0, 8.0) StateXY(9.0, 10.0)    StateXY(11.0, 2.0)]
# rget_state_df(state)

# Create a big matrix of StateXY objects
# np = 1000
# nt = 20000
# state = [StateXY(rand(), rand()) for _ in 1:np, _ in 1:nt]
# @time rget_state_df(state)