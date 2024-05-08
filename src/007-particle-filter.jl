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

A particle filtering algorithm that samples from f(X_t | {Y_1 ... Y_t}) for t ∈ 1:Tmax.

# Arguments:

- `timeline`: A Vector{DateTime} of time stamps that defines the time steps for the simulation. Time stamps must be equally spaced.
- `xinit`: A Vector{State} that defines the initial state(s) of the animal.
- `yobs`: A Dictionary of observations.
- `move`: A `ModelMove` instance.
- `n_move`: An integer that defines the number of attempts used to find a legal move. Particles are killed otherwise.
- `n_record`: An integer that defines the number of particles to record at each time step.
- `n_resample`: A number that defines the effective sample size at which to resample particles.

# Returns
A tuple with the following fields:
- `timeline`
- `state`
- `ess`
- `maxlp`
- `convergence`

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
