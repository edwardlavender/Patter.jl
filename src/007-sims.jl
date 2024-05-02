export sim_walk, sim_yobs

# Simulate movement path(s) 
# * xinit is a vector that defines the initial state(s)
# * The function simulates new movement path(s)
# * The output is a matrix as for particle_filter: 
# * rows: paths
# * columns: time steps
function sim_walk(; xinit = Vector, move::ModelMove, nt::Int64 = 1000)

    # Define output objects
    np = length(xinit)
    xout = Matrix{eltype(xinit)}(undef, np, nt);
    xout[:, 1] .= xinit
    
    # Run simulation
    for t in 2:nt
        for i in 1:np
            xout[i, t] = rmove(xout[i, t - 1], move, t, 1_000_000)[1]
        end
    end 

    xout

end

function initialise_dict(xinit::State, model::Vector, timestamp)
    tuple_type = [typeof(sim_obs(xinit, m, 1)) for m in model]
    tuple_model = [typeof(m) for m in model]
    tuples = [Tuple{tuple_type[i], tuple_model[i]} for i in 1:length(tuple_model)]
    Dict{eltype(timestamp), Vector{Union{tuples...}}}()
end

# Simulate observations
function sim_yobs(; paths::Matrix, models::Vector)
    
    #### Initialise a set of dictionaries (one per path)
    yobs_by_path = Dict()

    #### Iterate over paths & populate dictionaries 
    timesteps = 1:size(paths, 2)
    for i in 1:size(paths, 1)
        # Initialise dictionary 
        yobs_by_path[i] = []
        path = paths[i, :]
        yobs = initialise_dict(path[1], models, timesteps)
        # Simulate observations 
        for (t, state) in enumerate(path)
            yobs[timesteps[t]] = []
            for m in models
                push!(yobs[timesteps[t]], (sim_obs(state, m, t), m))
            end
        end
    
    push!(yobs_by_path[i], yobs)

    end 

    yobs_by_path

end


#### Define movement model
bathy = GeoArrays.read("/Users/lavended/Documents/work/projects/particle-filters/patter/projects/PatterDev.jl/data/bathy.tif")
move = ModelMoveXY(bathy, 
                   truncated(Gamma(1, 250.0), upper = 750.0), 
                   Uniform(-pi, pi))

#### Simulate a path
xinit = pf_init(StateXY[], move, 5, (707_000, 708_000), (6.267e6, 6.268e6));
paths = sim_walk(xinit = xinit, move = move, nt = 20000)

#### Simulate observations
models =  [ModelObsAcousticLogisStd(7.10347e5, 6.26696e6, 400.0, 100.0), 
          ModelObsAcousticLogisStd(7.08592e5, 6.26329e6, 400.0, 100.0), 
          ModelObsDepthUniform(bathy, 5.0, 5.0)]

obs = sim_yobs(paths = paths, models = models)

obs[1]