using Dates
using LogExpFunctions: logsumexp
using Base.Threads: @threads
using ProgressMeter: @showprogress

export pf_filter

# Particle filter
function pf_filter(
    ; timeline::Vector{DateTime}, 
    xinit::Vector,
    yobs::Dict,
    move::ModelMove, 
    n_move::Int = 1000,
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
                xnow[i], lwi = rmove(xpast[i], move, t, n_move)
                lw[i] += lwi
            end 
        end 

        #### Update (log) weights
        timestamp = timeline[t]
        if haskey(yobs, timestamp)
            @threads for (obs, model) in yobs[timestamp]
                for i in 1:np
                    if isfinite(lw[i])
                        lw[i] += log_prob_obs(xnow[i], model, t, obs)
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
        idx = pf_resample_lv(exp.(lw_norm), np)
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