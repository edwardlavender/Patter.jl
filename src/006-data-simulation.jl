export sim_path_walk

"""
# Simulate movement paths

This function simulates movement path(s) from a vector of initial states and movement model. 

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
function sim_path_walk(; xinit = Vector, move::ModelMove, nt::Int64 = 1000)

    # Define output objects
    np = length(xinit)
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
