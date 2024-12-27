using Random
using Dates
using LogExpFunctions: logsumexp
using Base.Threads: @threads
using ProgressMeter: @showprogress

export particle_filter, particle_filter_iter


"""
    Particles(states::Matrix, diagnostics::DataFrame, callstats::DataFrame)

# Fields

- `states`: A `Matrix` of [`State`](@ref)s:
    - Each row corresponds to a particle;
    - Each column corresponds to the `timestep`;
- `diagnostics`: A `DataFrame` of algorithm diagnostics:
    -   timestep: A `Vector{Int64}` of time steps;
    -   timestamp: A `Vector{DateTime}` of time stamps;
    -   ess: A `Vector{Float64}` that defines the effective sample size at each time step;
    -   maxlp:  `Vector{Float64}` that defines the maximum log weight at each time step (i.e., the maximum log-posterior, if resampling is implemented at every time step);
- `callstats`: A one-row `DataFrame` of call statistics:
    -   timestamp: A `DateTime` that define the start time of the function call;
    -   routine: A `String` that defines the algorithm;
    -   n_particle: An `Int` that defines the number of particles;
    -   n_iter: An `Int` or `NaN` that defines the number of iterations (trials);
    -   convergence: A `Boolian` that defines whether or not the algorithm reached the end of the `timeline`;
    -   time: A `Float64` that defines the duration (s) of the function call;

# See also

`Particles` objects are returned by:
    - [`particle_filter()`](@ref) 
    - [`two_filter_smoother()`](@ref)

"""
struct Particles
    states::Matrix{<:State}
    diagnostics::DataFrame
    callstats::DataFrame
end

# Create 'particles' structure in particle_filter() or two_filter_smoother()
function particulate(routine::String, 
                     timestamp::Dates.DateTime, 
                     timeline::Vector{Dates.DateTime},
                     states::Matrix{<:State},
                     ess::Vector{Float64}, 
                     maxlp::Vector{Float64}, 
                     n_particle::Int, 
                     n_iter::Union{Int, Float64}, # NaN in two_filter_smoother()
                     convergence::Bool)

    diagnostics = DataFrame(timestep  = collect(1:length(timeline)), 
                            timestamp = timeline, 
                            ess       = ess, 
                            maxlp     = maxlp)

    callstats = DataFrame(timestamp   = timestamp, 
                          routine     = routine,
                          n_particle  = n_particle,
                          n_iter      = n_iter,
                          convergence = convergence,
                          time        = call_duration(timestamp))
                          
    return Particles(states, diagnostics, callstats)

end 


"""
    resample(w::Vector{Float64}, n::Int)

(Internal) Given a weight vector `w`, resample a set of *indices* based on low-variance resampling algorithm from Thrun, Burgard and Fox's "Probabilistic Robotics".

# Arguments

- `w`: A `Vector{Float64}` of weights;
- `n`: An integer that defines the number of particles;

# Details

This is an internal function that implements systematic resampling in the particle filter (see [`particle_filter()`](@ref)) and smoothing algorithms (see [`two_filter_smoother()`](@ref)). Note that for large `n`, the function is not numerically stable.

# Returns

- An integer vector of indices;

# Example

```Julia
X = ["A", "B", "C", "D"]
w = [0, 0, 0.75, 0.25]

idx = resample(w, 12)
X[idx]
```

# Source

Code adapted from https://github.com/JuliaStats/StatsBase.jl/issues/124.

"""
function resample(w::Vector{Float64}, n::Int = length(w))
    if all(x -> x == 0 || isnan(x), w)
        error("All elements are zero/NaN.")
    end
    w = w ./ sum(w)
    idx = zeros(Int, n)
    r = rand() * 1/n
    c = w[1]
    i = 1
    for m in 1:n
        # Calculate the next sample point
        U = r + (m - 1) * (1 / n)
        # Find the first weight that puts us past the sample point
        while c < U && i < n
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

A particle filtering algorithm that samples from `f(X_t | {Y_1 ... Y_t}) for t ∈ 1:t`.

# Arguments (keywords)

- `timeline`: A `Vector{DateTime}` of ordered, regularly spaced time stamps that defines the time steps for the simulation;
- `xinit`: A `Vector` of [`State`](@ref) instances that defines the initial state(s) of the animal;
- `yobs`: A Dictionary of observations:
  - Dictionary keys should match elements in `timeline`;
  - Each element must be a `Vector` of `Tuple`s for that time step (one for each observation/sensor);
  - Each `Tuple` should contain (a) the observation and (b) the model parameters (that is, a [`ModelObs`](@ref) instance);
- `model_move`: A [`ModelMove`](@ref) instance:
    - The movement model describes movement from one time step to the next and therefore depends implicitly on the resolution of `timeline`;
    - The movement model should align with the [`State`](@ref) instances in `xinit`. For example, a two-dimensional state (`StateXY`) requires a corresponding movement model instance (i.e., `ModelMoveXY`);
- `n_move`: An integer that defines the number of attempts used to find a legal move;
    - All [`ModelMove`](@ref) sub-types contain a `map` field that defines the region(s) within which movements are allowed (see [`is_valid()`](@ref));
    - Each particle is moved up to `n_move` times, until a valid movement is simulated;
    - Particles that fail to generate a valid move are killed;
- `n_record`: An integer that defines the number of particles to record at each time step:
    - `n_record` particles are resampled at each time step and recorded in memory;
- `n_resample`: A number that defines the effective sample size for resampling:
    - Particles are resampled when the effective sample size <= `n_resample`;
- `t_resample`: `nothing`, an `integer` or a Vector of `integer`s that define the time step(s) at which to force resampling;
    - Particles are resampled at `t_resample` regardless of the effective sample size;
- `direction:` A `String` that defines the direction of the filter:
    - `"forward"` runs the filter forwards in time;
    - `"backward"` runs the filter backwards in time;

# Algorithm

## Initiation

The algorithm is initiated using a `Vector` of `n_particle` [`State`](@ref)s (`xinit`). See [`simulate_states_init()`](@ref) to simulate initial states for the filter.

## Movement

For every time step in the `timeline`, the internal function [`Patter.simulate_move()`](@ref) simulates the movement of particles away from previous [`State`](@ref)s into new [`State`](@ref)s using the movement model, as specified by `model_move`. [`Patter.simulate_move()`](@ref) is an iterative wrapper for a [`Patter.simulate_step()`](@ref) method that simulates a new [`State`](@ref) instance from the previous [`State`](@ref). [`Patter.simulate_move()`](@ref) implements [`Patter.simulate_step()`](@ref) iteratively until a legal move is found (or `n_move` is reached). For custom [`State`](@ref) or [`ModelObs`](@ref) sub-types, a corresponding [`Patter.simulate_step()`](@ref) method is required. Illegal moves are those that land in `NaN` locations on the `map` or, in the case of [`State`](@ref)s that include a depth (`z`) component, are below the depth of the seabed (see [`is_valid()`](@ref)). Particles that fail to generate legal moves are eventually killed by re-sampling (see below).

## Likelihood

Observations are used to weight simulated particles. To simulate observations for filtering, use [`simulate_yobs()`](@ref). To assemble real-world observations for filtering, see [`assemble_yobs()`](@ref). For each valid [`State`](@ref) and time stamp in `yobs`, the log-probability of each observation, given the [`State`](@ref), is evaluated via `Patter.logpdf_obs()`. For custom [`State`](@ref) or [`ModelObs`](@ref) sub-types, a corresponding `Patter.logpdf_obs()` method is required. The maximum weight across all particles (`maxlp`) is recorded at each time step as an algorithm diagnostic. (This metric can be intepreted as the maximum log-posterior if resampling is implemented at every time step.)

## Resampling

Particles are periodically re-sampled, with replacement, using the low-variance systematic re-sampling algorithm (via [`Patter.resample()`](@ref)), at time steps in `t_resample` or when the effective sample size is less than or equal to `n_resample`. This has the effect of eliminating impossible particles and duplicating likely ones.

The algorithm continues in this way, iterating over the `timeline`, simulating, weighting and (re)sampling particles. At each time step, `n_record` particles are saved in memory. If the function fails to converge, a warning is returned alongside the outputs up to that time step. Otherwise, the function will continue to the end of the time series.

## Multi-threading

The iteration over particles (i.e., simulated movements and likelihood evaluations) are multi-threaded.

## Convergence and diagnostics

Algorithm convergence is not guaranteed. The algorithm may reach a dead-end---a time step at which there are no valid locations into which the algorithm can step. This may be due to data errors, incorrect assumptions, insufficient sampling effort or poor tuning-parameter settings.

# Returns

- A [`Particles`](@ref) structure;

# See also

* [`State`](@ref), [`ModelMove`](@ref) and [`ModelObs`](@ref) for [`State`](@ref), movement model and observation model sub-types;
- [`simulate_yobs()`](@ref) and [`assemble_yobs()`](@ref) to prepare observations for the particle filter;
* [`Patter.simulate_step()`](@ref) and [`Patter.simulate_move()`](@ref) for the internal routines used to simulate new [`State`](@ref)s;
* `Patter.logpdf_obs()` methods to evaluate the log probability of observations;
* [`two_filter_smoother()`](@ref) to implement particle smoothing;

"""
function particle_filter(
    ; timeline::Vector{DateTime},
    xinit::Vector,
    yobs::Dict,
    model_move::ModelMove,
    n_move::Int = 100_000,
    n_record::Int = 1000,
    n_resample::Float64 = Float64(n_record),
    t_resample::Union{Nothing, Int, Vector{Int}} = nothing,
    direction::String = "forward")

    #### Define essential parameters
    # Call start
    call_start = now()
    # Number of time steps
    nt = length(timeline)
    # Number of particles
    np = length(xinit)
    # Number of recorded particles
    nr = n_record
    # Use t_resample
    do_t_resample = !isnothing(t_resample)

    #### Check user inputs
    # Check time line
    timeline = sort(timeline)
    check_timeline(timeline, collect(keys(yobs)))
    # Check particle numbers
    if nr > np
        error("The number of initial particles in `xinit` ($np) must be >= the number of recorded partices in `n_record` ($nr).")
    end

    #### Define filter direction
    if direction == "forward"
        start = 1; finish = nt;
        timesteps = start:finish
    elseif direction == "backward"
        start = nt; finish = 1;
        timesteps = start:-1:finish
    else
        error("`direction` must be \"forward\" or \"backward\".")
    end

    #### Define particle objects
    # Particles
    xpast = deepcopy(xinit)
    xnow  = deepcopy(xinit)
    xout  = Matrix{eltype(xinit)}(undef, nr, nt);
    # (log) weights
    lw = zeros(np)
    # Output particles
    xout = Matrix{eltype(xinit)}(undef, nr, nt);

    #### Define diagnostic objects
    # Output ESS vector
    ess = zeros(nt)
    # Output maxlp vector
    maxlp = zeros(nt)

    #### Run filter
    @showprogress desc = "Running filter..." for t in timesteps

        # println(t)

        #### Move particles & compute weights
        # * We iterate once over particles b/c this is thread safe
        timestamp            = timeline[t]
        has_obs_at_timestamp = haskey(yobs, timestamp)
        @threads for i in 1:np
            if isfinite(lw[i])
                # Move particles
                if t != start
                    xnow[i], lwi = simulate_move(xpast[i], model_move, t, n_move)
                    lw[i] += lwi
                end
                # Evaluate likelihoods
                if has_obs_at_timestamp && isfinite(lw[i])
                    for (obs, model) in yobs[timestamp]
                        lw[i] += logpdf_obs(xnow[i], model, t, obs)
                    end
                end
            end
        end

        #### Record diagnostics
        maxlp[t] = maximum(lw)

        #### Validate weights
        if !any(isfinite.(lw))
            # stop = ifelse(direction == "forward", t - 1, t + 1)
            stop = t
            pos  = sort([start, stop])
            pos  = pos[1]:pos[2]
            julia_warning("Weights from filter ($start -> $finish) are zero at time $t): returning outputs from $(minimum(pos)):$(maximum(pos)). Note that all (log) weights at $t are -Inf.")
            return particulate("filter: " * direction, call_start, 
                               timeline[pos], xout[:, pos], ess[pos], maxlp[pos], 
                               np, 1, false)
         end

        #### Resample particles
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
        do_resample = (do_t_resample && (t in t_resample)) || (ess[t] <= n_resample)
        if do_resample
            xpast .= xnow[idx]
            lw .= zero(Float64)
        else 
            xpast .= xnow
        end

    end

    return particulate("filter: " * direction, call_start, 
                       timeline, xout, ess, maxlp, 
                       np, 1, true)

end


"""
# Iterative particle filter

An interative implementation of the particle filter.

# Arguments (keywords)

- `...`: Keyword arguments passed to [`particle_filter`](@ref);
- `n_iter`: A integer that defines the maximum number of iterations (trials);

# Details

This function wraps [`particle_filter`](@ref). The filter is implemented up to `n_iter` times, or until convergence is achieved.

# Returns 

- A [`Particles`](@ref) structure;

# See also

* [`particle_filter`](@ref) implements the particle filter;

"""
function particle_filter_iter(
    ; timeline::Vector{DateTime},
    xinit::Vector,
    yobs::Dict,
    model_move::ModelMove,
    n_move::Int = 100_000,
    n_record::Int = 1000,
    n_resample::Float64 = Float64(n_record),
    t_resample::Union{Nothing, Int, Vector{Int}},
    n_iter::Int64 = 1,
    direction::String = "forward")

    # Run filter iteratively 
    call_start = now()
    iter       = 0
    n_iter     = max(1, n_iter)
    out        = nothing
    while iter < n_iter
        iter = iter + 1
        out  = particle_filter(timeline   = timeline,
                               xinit      = xinit,
                               yobs       = yobs,
                               model_move = model_move,
                               n_move     = n_move,
                               n_record   = n_record,
                               n_resample = n_resample,
                               t_resample = t_resample,
                               direction  = direction)
        if out.callstats.convergence[1]
            break
        end 
    end
    
    # Define outputs
    out.callstats.timestamp .= call_start
    out.callstats.n_iter    .= iter
    out.callstats.time      .= call_duration(call_start)
    return out

end