using Base.Threads: @threads
using Dates
using ProgressMeter: @showprogress
using JLD2
using LogExpFunctions: logsumexp
using Random

export particle_filter


"""
    Particles(states::Union{Nothing, Matrix{<:State}}, diagnostics::DataFrame, callstats::DataFrame)

# Fields

- `states`: (optional) A `Matrix` of [`State`](@ref)s:
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
    -   convergence: A `Boolian` that defines convergence;
    -   time: A `Float64` that defines the duration (s) of the function call;

# Details

* `states` is `nothing` if [`particle_filter()`](@ref) or [`particle_smoother_two_filter()`](@ref) are implemented with `batch`ing. 
* `convergence` is defined as follows:
    - In [`particle_filter()`](@ref), `convergence` defines whether or not the filter reached the end of the `timeline`;
    - In [`particle_smoother_two_filter()`](@ref), `convergence` defines whether or not correct smoothing was achieved on at least 95 % of time steps. 'Correct smoothing' is possible when there at at least some valid moves between the subset of recorded particles on the backward filter and those on the forward filter (for the previous time step). 

# See also

`Particles` objects are returned by:
- [`particle_filter()`](@ref)
- [`particle_smoother_two_filter()`](@ref)

"""
struct Particles
    states::Union{Nothing, Matrix{<:State}}
    diagnostics::DataFrame
    callstats::DataFrame
end

# Create 'particles' structure in particle_filter() or particle_smoother_two_filter()
function particulate(routine::String, 
                     timestamp::Dates.DateTime, 
                     timeline::Vector{Dates.DateTime},
                     states::Union{Nothing, Matrix{<:State}},
                     ess::Vector{Float64}, 
                     maxlp::Vector{Float64}, 
                     n_particle::Int, 
                     n_iter::Union{Int, Float64}, # NaN in particle_smoother_two_filter()
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

This is an internal function that implements systematic resampling in the particle filter (see [`particle_filter()`](@ref)) and smoothing algorithms (see [`particle_smoother_two_filter()`](@ref)). Note that for large `n`, the function is not numerically stable.

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


# Internal particle filter function 
function _particle_filter(
    ; timeline::Vector{DateTime},
    xinit::Vector{<:State},
    yobs::Dict,
    model_move::ModelMove,
    n_move::Int = 100_000,
    n_record::Int = 1000,
    n_resample::Float64 = Float64(n_record),
    t_resample::Union{Nothing, Int, Vector{Int}} = nothing,
    direction::String = "forward", 
    batch::Union{Nothing, Vector{String}} = nothing)

    #### Define essential parameters
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

    #### Handle batches
    do_batch = !isnothing(batch)
    nb       = do_batch ? length(batch) : 1

    #### Define filter direction
    if direction == "forward"
        start = 1; finish = nt;
        timesteps = start:finish
        if do_batch
            batch = sort(batch)
        end 
    elseif direction == "backward"
        start = nt; finish = 1;
        timesteps = start:-1:finish
        if do_batch
            batch = sort(batch, rev = true)
        end 
    else
        error("`direction` must be \"forward\" or \"backward\".")
    end
    # Batch time steps
    timesteps = collect(timesteps)
    timesteps_by_batch = split_indices(timesteps, nb)

    #### Define filter
    # Particles
    xpast = deepcopy(xinit)
    xnow  = deepcopy(xinit)
    xout  = Matrix{eltype(xinit)}(undef, nr, length(timesteps_by_batch[1])) 
    # (log) weights
    lw = zeros(np)
    # Output ESS vector
    ess = fill(NaN, nt)
    # Output maxlp vector
    maxlp = fill(NaN, nt)

    #### Run filter
    for b in 1:nb

        # Define output particles object
        xout = Matrix{eltype(xinit)}(undef, nr, length(timesteps_by_batch[b]))

        # Define time steps and indicies for iteration 
        timesteps_for_batch = timesteps_by_batch[b]
        indices_for_batch   = collect(1:length(timesteps_for_batch))
        if direction == "backward"
            indices_for_batch = reverse(indices_for_batch)
        end 

        # Run filter
        @showprogress desc = "Running filter..." for (i, t) in zip(indices_for_batch, timesteps_for_batch)

            # println(t)

            #### Move particles & compute weights
            # * We iterate once over particles b/c this is thread safe
            timestamp            = timeline[t]
            has_obs_at_timestamp = haskey(yobs, timestamp)
            @threads for j in 1:np
                if isfinite(lw[j])
                    # Move particles
                    if t != start
                        xnow[j], lwi = simulate_move(xpast[j], model_move, t, n_move)
                        lw[j] += lwi
                    end
                    # Evaluate likelihoods
                    if has_obs_at_timestamp && isfinite(lw[j])
                        for (obs, model) in yobs[timestamp]
                            lw[j] += logpdf_obs(xnow[j], model, t, obs)
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
                if do_batch
                    if direction == "forward"
                        @save batch[b] xfwd = xout
                    else
                        @save batch[b] xbwd = xout
                    end 
                    xout = nothing
                end 
                return (states = xout, ess = ess, maxlp = maxlp, convergence = false)
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
            xout[:, i] .= xnow[idx[1:nr]]
            # Optionally resample particles for next iteration
            do_resample = (do_t_resample && (t in t_resample)) || (ess[t] <= n_resample)
            if do_resample
                xpast .= xnow[idx]
                lw .= zero(Float64)
            else 
                xpast .= xnow
            end

        end

        # Write outputs to batch for file
        if do_batch
            if direction == "forward"
                @save batch[b] xfwd = xout
            else
                @save batch[b] xbwd = xout
            end 
            xout = nothing
        end 

    end 

    return (states = xout, ess = ess, maxlp = maxlp, convergence = true)


end

"""
    particle_filter(; timeline::Vector{DateTime},
                      xinit::Vector{<:State},
                      yobs::Dict,
                      model_move::ModelMove,
                      n_move::Int = 100_000,
                      n_record::Int = 1000,
                      n_resample::Float64 = Float64(n_record),
                      t_resample::Union{Nothing, Int, Vector{Int}},
                      n_iter::Int64 = 1,
                      direction::String = "forward", 
                      batch::Union{Nothing, Vector{String}} = nothing)

A particle filtering algorithm that samples from `f(s_t | y_{1:t})` for `t ∈ 1:t`.

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
- `n_iter`: A integer that defines the maximum number of iterations (trials);
- `direction:` A `String` that defines the direction of the filter:
    - `"forward"` runs the filter forwards in time;
    - `"backward"` runs the filter backwards in time;
- (optional) `batch`: A Vector of `.jld2` file paths for particles (see Memory Management);

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

## Memory management 

By default, `n_record` particles at each time step are retained in memory. If `batch` is provided, the `timeline` is split into `length(batch)` batches. The filter still moves along the whole `timeline`, but only records the particles for the current batch in memory. At the end of each batch, the particles for that batch are written to file. This reduces total memory demand. 

`batch` file paths are sorted alphanumerically if `direction = "forward"` and in reverse order if `direction = "backward"`. For example: 
* If you have a `timeline` of 10 time steps, `direction = "forward"` and `batch = ["fwd-1.jld2", "fwd-2.jld2", "fwd-3.jld2"]`, `fwd-1.jld2`, `fwd-2.jld2` and `fwd-3.jld2` contain the particle matrices for time steps `[1, 2, 3]`, `[4, 5, 6]` and [`7, 8, 9, 10`], respectively. 
* If you have a `timeline` of 10 time steps, `direction = "backward"` and `batch = ["bwd-1.jld2", "bwd-2.jld2", "bwd-3.jld2"]`, `bwd-1.jld2`, `bwd-2.jld2` and `bwd-3.jld2` similarly contain the particle matrices for time steps `[1, 2, 3]`, `[4, 5, 6]` and [`7, 8, 9, 10`], respectively. 

## Convergence and diagnostics

Algorithm convergence is not guaranteed. The algorithm may reach a dead-end---a time step at which there are no valid locations into which the algorithm can step. This may be due to data errors, incorrect assumptions, insufficient sampling effort or poor tuning-parameter settings.

# Returns

- A [`Particles`](@ref) structure;

# See also

* [`State`](@ref), [`ModelMove`](@ref) and [`ModelObs`](@ref) for [`State`](@ref), movement model and observation model sub-types;
- [`simulate_yobs()`](@ref) and [`assemble_yobs()`](@ref) to prepare observations for the particle filter;
* [`Patter.simulate_step()`](@ref) and [`Patter.simulate_move()`](@ref) for the internal routines used to simulate new [`State`](@ref)s;
* `Patter.logpdf_obs()` methods to evaluate the log probability of observations;
* [`particle_smoother_two_filter()`](@ref) to implement particle smoothing;

"""
function particle_filter(
    ; timeline::Vector{DateTime},
    xinit::Vector{<:State},
    yobs::Dict,
    model_move::ModelMove,
    n_move::Int = 100_000,
    n_record::Int = 1000,
    n_resample::Float64 = Float64(n_record),
    t_resample::Union{Nothing, Int, Vector{Int}},
    n_iter::Int64 = 1,
    direction::String = "forward", 
    batch::Union{Nothing, Vector{String}} = nothing)

    # Run filter iteratively 
    call_start = now()
    iter       = 0
    n_iter     = max(1, n_iter)
    out        = nothing
    while iter < n_iter
        iter = iter + 1
        out  = _particle_filter(timeline   = timeline,
                                xinit      = xinit,
                                yobs       = yobs,
                                model_move = model_move,
                                n_move     = n_move,
                                n_record   = n_record,
                                n_resample = n_resample,
                                t_resample = t_resample,
                                direction  = direction, 
                                batch      = batch)
        if out.convergence[1]
            break
        end 
    end

    return particulate("filter: " * direction, call_start, 
                       timeline, out.states, out.ess, out.maxlp, 
                       length(xinit), iter, out.convergence)

end