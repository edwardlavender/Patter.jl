using Random
using Dates
using LogExpFunctions: logsumexp
using Base.Threads: @threads
using ProgressMeter: @showprogress

export particle_filter


# Check user inputs

function check_timeline_entries(t_sim, t_obs)
    issubset(t_obs, t_sim) || error("There are time stamps in `yobs` not found in `timeline`.")
    nothing
end 

function check_timeline_spacing(t_sim)
    all(diff(t_sim) .== first(diff(t_sim))) || error("`timeline` must be a sequence of requally spaced time steps.")
    nothing
end 

function check_timeline(t_sim, t_obs)
    check_timeline_entries(t_sim, t_obs)
    check_timeline_spacing(t_sim)
    nothing
end 


"""
    # Systematic resampling algorithm

    Given the weight vector w, resample a set of *indices* based on low-variance resampling algorithm from Thrun, Burgard, and Fox's "Probabilistic Robotics".

    # Source

    Code adapted from https://github.com/JuliaStats/StatsBase.jl/issues/124.

    # Example

    ```Julia
    X = ["A", "B", "C", "D"]
    w = [0, 0, 0.75, 0.25]

    idx = resample(w, 12)
    X[idx]
    ```
    """
function resample(w::Vector{Float64}, n::Int = length(w))
    w = w ./ sum(w)
    idx = zeros(Int, n)
    r = rand() * 1/n
    c = w[1]
    i = 1
    for m in 1:n
        # Calculate the next sample point
        U = r + (m - 1) * (1 / n)
        # Find the first weight that puts us past the sample point
        while c < U
            i += 1
            c = c + w[i]
        end
        idx[m] = i
    end
    Random.shuffle!(idx)
    idx
end


"""
# Particle filter

A particle filtering algorithm that samples from `f(X_t | {Y_1 ... Y_t}) for t ∈ 1:Tmax`.

# Arguments

- `timeline`: A Vector{DateTime} of ordered, regularly spaced time stamps that defines the time steps for the simulation.
- `xinit`: A Vector of [`State`](@ref) instances that defines the initial state(s) of the animal.
- `yobs`: A Dictionary of observations (see [`assemble_yobs()`](@ref)):  
        - Dictionary keys should match elements in `timeline`. 
        - Each element must be a `Vector` of `Tuple`s for that time step (one for each observation/sensor).  
        - Each `Tuple` should contain (a) the observation and (b) the model parameters (that is, a `ModelObs` instance).
- `move`: A [`ModelMove`](@ref) instance.
    - The movement model describes movement from one time step to the next and therefore depends implicitly on the resolution of `timeline`.
    - The movement model should align with the [`State`] instances in `.xinit`. For example, a 2-dimensional state ([`StateXY`](@ref)) requires a corresponding movement model instance (i.e., [`ModelMoveXY`](@ref)). 
- `n_move`: An integer that defines the number of attempts used to find a legal move. 
    - All [`ModelMove`](@ref) sub-types contain a `map` field that defines the region(s) within which movements are allowed.
    - Each particle is moved up to `n_move` times, until a valid movement is simulated. 
    - Particles that fail to generate a valid move are killed. 
- `n_record`: An integer that defines the number of particles to record at each time step.
    - `n_record` particles are resampled at each time step and recorded in memory. 
- `n_resample`: A number that defines the effective sample size for resampling.
    - Particles are resampled when the effective sample size <= `n_resample`.
- `direction:` A `String` that defines the direction of the filter.
    - "forward" runs the filter forwards in time;
    - "backward" runs the filter backwards in time;

# Algorithm

## Initiation 
The algorithm is initiated using a `Vector` of `.n_particle` [`State`](@ref)s (`.xinit`). 

## Movement 

For every time step in the `timeline`, the internal function `Patter.simulate_move()` simulates the movement of particles away from previous states into new states using the movement model specified by `.model_move`. `Patter.simulate_move()` is an iterative wrapper for a `Patter.simulate_step()` method that simulates a new [`State`](@ref) instance from the previous [`State`](@ref), using the movement model. `Patter.simulate_move()` implements `Patter.simulate_step()` iteratively until a legal move is found (or `.n_move` is reached). Illegal moves are those that land in `NaN` locations on the `map` or, in the case of states that include a depth (`z`) component, are below the depth of the seabed. Particles that fail to generate legal moves are eventually killed by re-sampling (see below).
    
## Likelihood 

For each valid [`State`](@ref) and time stamp in `yobs`, the log-probability of each observation, given the [`State`](@ref), is evaluated via `Patter.logpdf_obs()`. The maximum log-probability across all particles is recorded at each time step as an algorithm diagnostic.

## Resampling 

Particles are periodically re-sampled, with replacement, using the low-variance systematic re-sampling algorithm (via [`resample()`](@ref)), when the effective sample size is less than or equal to `.n_resample`. This has the effect of eliminating impossible particles and duplicating likely ones.

The algorithm continues in this way, iterating over the `timeline`, simulating, weighting and (re)sampling particles. At each time step, `.n_record` particles are saved in memory. If the function fails to converge, a [`warning`] is returned alongside the outputs up to that time step. Otherwise, the function will continue to the end of the time series.

# Multi-threading

The iteration over particles (i.e., simulated movements and likelihood evaluations) are multi-threaded. 

# Convergence and diagnostics

Algorithm convergence is not guaranteed. The algorithm may reach a dead-end---a time step at which there are no valid locations into which the algorithm can step. This may be due to data errors, incorrect assumptions, insufficient sampling effort or poor tuning-parameter settings.

# Returns
A `NamedTuple` with the following fields:
- `timesteps`: An `Vector{Int64}` of time steps;
- `timestamps`: The `timeline`;
- `direction`: The `direction`;
- `state`: A `Matrix` of [`State`](@ref)s:
    - Each row corresponds to a particle; 
    - Each column corresponds to the `timestep`;
- `ess`: A `Vector{Float64}` that defines the effective sample size at each time step;
- `maxlp`: A `Vector{Float64}` that defines the maximum log-posterior at each time step;
- `convergence`: A `Bool` that defines whether or not the algorithm reached the end of the `timeline`;
"""
function particle_filter(
    ; timeline::Vector{DateTime},
    xinit::Vector,
    yobs::Dict,
    move::ModelMove,
    n_move::Int = 100_000,
    n_record::Int = 1000,
    n_resample::Float64 = n_record * 0.5, 
    direction::String = "forward")

    #### Define essential parameters
    # Number of time steps
    nt = length(timeline)
    # Number of particles
    np = length(xinit)
    # Number of recorded particles
    nr = n_record

    #### Check user inputs
    # Check time line
    timeline = sort(timeline)
    check_timeline(timeline, keys(yobs))
    # Check particle numbers
    if nr > np
        error("The number of initial particles in `xinit` ($np) must be >= the number of recorded partices in `n_record` ($nr).")
    end

    #### Define filter direction
    if direction == "forward"
        start = 1
        timesteps = 2:nt
    elseif direction == "backward"
        start = nt
        timesteps = (nt - 1):-1:1
    else 
        error("`direction` must be \"forward\" or \"backward\".")
    end 

    #### Define particle objects
    # Particles
    xpast = deepcopy(xinit)
    xnow = deepcopy(xinit)
    xout = Matrix{eltype(xinit)}(undef, nr, nt);
    # (log) weights
    lw = zeros(np)
    # Output particles
    xout = Matrix{eltype(xinit)}(undef, nr, nt);
    xout[:, start] .= xinit[1:nr]

    #### Define diagnostic objects
    # Output ESS vector
    ess = zeros(nt)
    ess[start] = np
    # Output maxlp vector
    maxlp = zeros(nt)
    maxlp[start] = NaN
   
    #### Run filter
    @showprogress desc = "Running filter..." for t in timesteps

        # println(t)

        #### Move particles & compute weights
        # * We iterate once over particles b/c this is thread safe
        timestamp            = timeline[t]
        has_obs_at_timestamp = haskey(yobs, timestamp)
        yobs_at_timestamp    = yobs[timestamp]
        @threads for i in 1:np
            if isfinite(lw[i])
                # Move particles
                xnow[i], lwi = simulate_move(xpast[i], move, t, n_move)
                lw[i] += lwi
                # Evaluate likelihoods
                if has_obs_at_timestamp && isfinite(lw[i])
                    for (obs, model) in yobs_at_timestamp
                        lw[i] += logpdf_obs(xnow[i], model, t, obs)
                    end
                end 
            end
        end

        #### Record diagnostics
        maxlp[t] = maximum(lw)

        #### (optional) Resample particles
        # Validate weights
        if !any(isfinite.(lw))
           @warn "Weights from filter ($start -> $t) are zero at time $t: returning outputs up to $(t - 1)."
           pos = sort([start, t - 1])
           pos = pos[1]:pos[2]
           return (timesteps    = collect(pos), 
                   timestamps   = timeline[pos], 
                   state        = xout[:, pos], 
                   direction    = direction, 
                   ess          = ess[pos], 
                   maxlp        = maxlp[pos], 
                   convergence  = false)
        end
        # Normalise weights
        lw_norm = lw .- logsumexp(lw)
        # Evaluate ESS
        ess[t] = exp(-logsumexp(2 * lw_norm))
        # Compute resampling indices
        idx = resample(exp.(lw_norm), np)
        # Record (subset) of resampled particles
        # (A deep copy is implicitly made here via the subsetting)
        xout[:, t] .= xnow[idx[1:nr]]
        # Optionally resample particles for next iteration
        if ess[t] <= n_resample
            xpast .= xnow[idx]
            lw .= zero(Float64)
        else
            xpast .= xnow
        end

    end

    (
        timesteps   = collect(1:nt), 
        timestamps  = timeline, 
        state       = xout, 
        direction   = direction, 
        ess         = ess, 
        maxlp       = maxlp, 
        convergence = true
    )

end