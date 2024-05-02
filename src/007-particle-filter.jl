using Random
using Dates
using LogExpFunctions: logsumexp
using Base.Threads: @threads
using ProgressMeter: @showprogress

export particle_filter

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
        U = r + (m - 1) * (1/n)
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
- `resample_ess`: A number that defines the effective sample size at which to resampler particles.

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
    resample_ess::Int = n_record * 0.5)

    #### Define essential parameters
    # Number of time steps
    nt = length(timeline)
    # Number of particles
    np = length(xinit)
    # Number of recorded particles
    nr = n_record

    #### Check user inputs
    if nr > np
        error("The number of initial particles in `xinit` ($np) must be >= the number of recorded partices in `n_record` ($nr).")
    end

    #### Initiate filter
    # Particles
    xpast = deepcopy(xinit)
    xnow = deepcopy(xinit)
    # (log) weights
    lw = zeros(np)
    # Output ESS vector
    ess = zeros(nt)
    ess[1] = np
    # Output maxlp vector
    maxlp = zeros(nt)
    maxlp[1] = NaN
    # Output particles
    xout = Matrix{eltype(xinit)}(undef, nr, nt);
    xout[:, 1] .= xinit[1:nr]

    #### Select direction
    # TO DO

    #### Run filter
    @showprogress desc = "Computing..." for t in 2:nt

        # println(t)

        #### Move particles
        # * We thread over particles & separately over likelihoods
        # * This appears to be faster than threading once over particles
        timestamp = timeline[t]
        @threads for i in 1:np
            if isfinite(lw[i])
                xnow[i], lwi = simulate_move(xpast[i], move, t, n_move)
                lw[i] += lwi
            end
        end

        #### Update (log) weights
        if haskey(yobs, timestamp)
            @threads for (obs, model) in yobs[timestamp]
                for i in 1:np
                    if isfinite(lw[i])
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
           @warn "All weights are zero at time $t: returning outputs up to this time."
           return (timestamp = timeline[1:t], state = xout[:, 1:t], ess = ess[1:t], maxlp = maxlp[1:t], convergence = false)
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
        if ess[t] <= resample_ess
            xpast .= xnow[idx]
            lw .= zero(Float64)
        else
            xpast .= xnow
        end

    end

    (timestamp = timeline, state = xout, ess = ess, maxlp = maxlp, convergence = true)

end
