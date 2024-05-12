using Base.Threads: @threads
using ProgressMeter: @showprogress

export two_filter_smoother

"""
Two filter smoother

Fearnhead, P., Wyncoll, D., Tawn, J., 2010. A sequential smoothing
algorithm with linear computational cost. Biometrika 97,
447–464. https://doi.org/10.1093/biomet/asq013

"""
function two_filter_smoother(;xfwd::Matrix, xbwd::Matrix, move::ModelMove, box = nothing, nMC::Int = 100)

    #### Check inputs
    size(xfwd) == size(xbwd) || error("Forward and backward sample do not match!")

    #### Set up
    # Dimension of input state
    zdim = hasfield(typeof(xfwd[1]), :z)
    # Matrix for smoothed particles
    xsmooth = similar(xfwd)
    xsmooth[:, 1] .= xbwd[:, 1]
    xsmooth[:, end] .= xfwd[:, end]
    np, nt = size(xsmooth)
    # Vectors for weights and ESS
    w = ones(np)
    ess = zeros(nt)
    # Least recently used cache 
    cache = LRU{eltype(xfwd), Float64}(maxsize = np)

    #### Run two-filter formula
    @showprogress desc = "Running two-filter smoother..." for t in 2:(nt - 1)
        @threads for k in 1:np
            for j in 1:np
                # Evaluate probability density of movement between locations (i.e., the weight)
                w[k] += exp(logpdf_move(xbwd[k, t], xfwd[j, t - 1], zdim, move, t, box, nMC, cache))
            end
        end
        # Normalise weights & evaluate ESS
        w .= w ./ sum(w)
        ess[t] = 1 / sum(abs2, w)
        # Resample particles & store
        idx = resample(w, np)
        xsmooth[:, t] .=  xbwd[idx, t]
    end

    (state = xsmooth, ess = ess)

end