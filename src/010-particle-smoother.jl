using Base.Threads: @threads
using ProgressMeter: @showprogress

export two_filter_smoother
export backward_smoother

"""
    two_filter_smoother(; timeline::Vector{DateTime}, xfwd::Matrix, xbwd::Matrix, model_move::ModelMove, vmap, n_sim::Int)

A two-filter particle smoother that samples from `f(X_t | {Y_1 ... Y_T}) for t ∈ 1:T`.

# Arguments (keywords)

- `timeline`: A `Vector{DateTime}` of ordered, regularly spaced time stamps that defines the time steps for the simulation;
- `xfwd`: A `Matrix` of [`State`](@ref)s from the forward filter (see [`particle_filter()`](@ref));
- `xbwd`: A `Matrix` of [`State`](@ref)s from the backward filter (see [`particle_filter()`](@ref));
- `model_move`: A [`ModelMove`](@ref) instance;
- `vmap`: (optional) A `GeoArray` that defines the 'validity map' (see [`logpdf_move()`](@ref));
- `n_sim`: An integer that defines the number of Monte Carlo simulations (see [`logpdf_move()`](@ref));

# Details

[`two_filter_smoother()`](@ref) smooths particles from the particle filter (see [`particle_filter()`](@ref)). The `timeline` from the particle filter should be supplied as well as a `Matrix` of particles from a forward run and a backward run. The two filter smoother works by iteratively resampling particles in line with the probability density of movement between particles from the backward filter at time `t` and particles from the forward filter at time `t - 1`. [`logpdf_move()`](@ref) is an internal function that evaluates the log probability of a movement step between particles. This function wraps the [`logpdf_step()`](@ref) generic. Methods are provided for built-in [`State`](@ref) and [`ModelMove`](@ref) sub-types. To use custom sub-types, a corresponding [`logpdf_step()`](@ref) method should be provided. In [`two_filter_smoother()`](@ref), the `vmap` and `n_sim` arguments support the calculate of probability densities (see [`logpdf_move()`](@ref)). 

# Returns 

- A `NamedTuple`, of the same format as returned by [`particle_filter()`](@ref), with the following fields:
    - `timesteps`
    - `timestamps`
    - `direction`: `nothing`
    - `state`
    - `ess`
    - `maxlp`: `NaN`
    - `convergence`: `true`
    - `trial`: NaN

# See also

- [`particle_filter()`](@ref) implements the particle filter;
- [`logpdf_step()`](@ref), [`logpdf_move_normalisation()`](@ref) and [`logpdf_move()`](@ref) evaluate the log probability (density) of movement between two [`State`](@ref)s;
- [`two_filter_smoother()`](@ref) implements the two-filter particle smoother;

# Source 

Fearnhead, P., Wyncoll, D., Tawn, J., [2010](https://doi.org/10.1093/biomet/asq013). A sequential smoothing algorithm with linear computational cost. Biometrika 97, 447–464.

"""
function two_filter_smoother(;timeline::Vector{DateTime}, xfwd::Matrix, xbwd::Matrix, model_move::ModelMove, vmap = nothing, n_sim::Int = 100)

    #### Check inputs
    size(xfwd) == size(xbwd) || error("Forward and backward sample do not match!")

    #### Set up
    # Identify dimension of input state
    zdim = hasfield(typeof(xfwd[1]), :z)
    # Use vmap, if supplied
    # * vmap may be supplied for 'horizontal' (2D) movement models 
    # * If vmap is supplied, we update we update xfwd.map_value to 0.0 or 1.0 
    # * If we are within mobility of the coastline, map_value -> 0.0
    # * map_value = 0.0 indicates that not all moves from that point are valid 
    # * This is a flag that we need to simulate the normalisation constant in logpdf_move()
    if !isnothing(vmap)
        xbwd = edit_map_value(xbwd, vmap)
    end 
    # Define smoothed particles matrix 
    # (rows: particles; columns: time steps)
    xout = similar(xfwd)
    xout[:, 1] .= xbwd[:, 1]
    xout[:, end] .= xfwd[:, end]
    # Define particle indices
    # * If the two filters are incompatible (all weights zero), 
    # * ... we sample n_fwd_half and n_bwd_half particles from  
    # * ... the forward and backward filters respectively
    np, nt = size(xout)
    half = (np ÷ 2)
    n_fwd_half = length(1:half)
    n_bwd_half = length((half + 1):np)
    warn_zero_weights = false 
    # Initialise ESS vector
    ess     = zeros(nt)
    ess[1]  = np # = 1 / sum(abs2, w) given equal weights from filter
    ess[nt] = np 
    # Set up LRU cache 
    # * This caches movement normalisation constants for s_{i, t} in logpdf_move()
    cache = LRU{eltype(xfwd), Float64}(maxsize = np)

    #### Run smoothing
    @showprogress desc = "Running two-filter smoother..." for t in 2:(nt - 1)

        # Compute weights  
        w = zeros(np)
        @threads for k in 1:np
            for j in 1:np
                # Evaluate probability density of movement between locations (i.e., the weight)
                w[k] += exp(logpdf_move(xbwd[k, t], xfwd[j, t - 1], zdim, model_move, t, vmap, n_sim, cache))
            end
        end
        
        # Validate weights & implement resampling 
        if any(w .> 0)
            # (A) If there are positive weights, normalise, compute ESS & resample
            w .= w ./ sum(w)
            ess[t] = 1 / sum(abs2, w)
            # Resample particles from xbwd & store
            idx = resample(w, np)
            xout[:, t] .=  xbwd[idx, t]
        else 
            # (B) If all weights are zero, set a warning flag 
            # * For speed, the warning is thrown only once outside the for loop 
            # * (as it may occur at multiple time steps)
            warn_zero_weights = true
            # Set ESS[t] = NaN
            ess[t] = NaN
            # Sample 50 % of particles from forward filter
            xout[1:half , t] = xfwd[rand(1:np, n_fwd_half), t]
            # Sample 50 % of particles from backward filter 
            xout[(half + 1):end, t] = xbwd[rand(1:np, n_bwd_half), t] 

        end 

    end

    #### Implement post-processing 
    # Warn for smoothing failures
    if warn_zero_weights
        nan_count = count(isnan, ess)
        nan_perc  = round(nan_count / length(ess) * 100, digits = 2)
        julia_warning("All smoothing weights (from xbwd[k, t] to xfwd[j, t - 1]) are zero at $nan_count time step(s) ($nan_perc %).")
    end 
    # Reset map_value
    xout = edit_map_value(xout, model_move.map)

    #### Return outputs
    # Follow particle_filter() format
    (
        timesteps   = collect(1:nt), 
        timestamps  = timeline, 
        state       = xout, 
        direction   = nothing, 
        ess         = ess, 
        maxlp       = NaN, 
        convergence = true, 
        trials      = NaN
    )

end

"""
    backward_smoother(;timeline::Vector{DateTime}, xfwd::Matrix, model_move::ModelMove, vmap, n_sim::Int)

The backward particle smoother that samples from `f(X_t | {Y_1 ... Y_T}) for t ∈ 1:T`.

# Arguments (keywords)

- `timeline`: A `Vector{DateTime}` of ordered, regularly spaced time stamps that defines the time steps for the simulation;
- `xfwd`: A `Matrix` of [`State`](@ref)s from the forward filter (see [`particle_filter()`](@ref));
- `model_move`: A [`ModelMove`](@ref) instance;
- `vmap`: (optional) A `GeoArray` that defines the 'validity map' (see [`logpdf_move()`](@ref));
- `n_sim`: An integer that defines the number of Monte Carlo simulations (see [`logpdf_move()`](@ref));

# Details

[`backward_smoother()`](@ref) smooths particles from the particle filter (see [`particle_filter()`](@ref)). The `timeline` from the particle filter should be supplied as well as a `Matrix` of particles from a forward run. The backward smoother works by iteratively resampling particles in line with the probability density of movement between the particles at time `t` and time `t+1`. [`logpdf_move()`](@ref) is an internal function that evaluates the log probability of a movement step between particles. This function wraps the [`logpdf_step()`](@ref) generic. Methods are provided for built-in [`State`](@ref) and [`ModelMove`](@ref) sub-types. To use custom sub-types, a corresponding [`logpdf_step()`](@ref) method should be provided. In [`backward_smoother`](@ref), the `vmap` and `n_sim` arguments support the calculate of probability densities (see [`logpdf_move()`](@ref)). 

# Returns 

- A `NamedTuple`, of the same format as returned by [`particle_filter()`](@ref), with the following fields:
    - `timesteps`
    - `timestamps`
    - `direction`: `nothing`
    - `state`
    - `ess`
    - `maxlp`: `NaN`
    - `convergence`: `true`
    - `trial`: NaN

# See also

- [`particle_filter()`](@ref) implements the particle filter;
- [`logpdf_step()`](@ref), [`logpdf_move_normalisation()`](@ref) and [`logpdf_move()`](@ref) evaluate the log probability (density) of movement between two [`State`](@ref)s;
- [`two_filter_smoother()`](@ref) implements the two-filter particle smoother;

"""
function backward_smoother(;timeline::Vector{DateTime}, xfwd::Matrix, model_move::ModelMove, vmap = nothing, n_sim::Int = 100)

    #### Set up
    # Identify dimension of input state
    zdim = hasfield(typeof(xfwd[1]), :z)
    # Use vmap, if supplied
    # * vmap may be supplied for 'horizontal' (2D) movement models 
    # * If vmap is supplied, we update we update xfwd.map_value to 0.0 or 1.0 
    # * If we are within mobility of the coastline, map_value -> 0.0
    # * map_value = 0.0 indicates that not all moves from that point are valid 
    # * This is a flag that we need to simulate the normalisation constant in logpdf_move()
    if !isnothing(vmap)
        xfwd = edit_map_value(xfwd, vmap)
    end 
    # Define smoothed particles matrix 
    # (rows: particles; columns: time steps)
    xout = similar(xfwd)
    xout[:, end] .= xfwd[:, end]
    np, nt = size(xout)
    # Initialise weights
    # * Particles from the forward filter have uniform weights (w) thanks to resampling
    # * At the last (first smoothing) time step, smoothed weights (ws) are also uniform 
    lw = log.(ones(np) ./ np)            
    ws = ones(np) ./ np
    # Initialise ESS vector
    ess = zeros(Float64, nt)
    ess[nt] = np # = 1 / sum(abs2, ws) 
    # Set up LRU cache 
    # * This caches movement normalisation constants for s_{i, t} in logpdf_move()
    cache = LRU{eltype(xfwd), Float64}(maxsize = np)

    #### Run smoothing
    @showprogress desc = "Running backward smoother..." for t in (nt-1):-1:1
   
        # Compute pairwise densities f(s_{j, t + 1} | f(s_{i, t})  * w
        # * We use a matrix as densities may be directional 
        w_i_to_j = zeros(Float64, np, np)
        @threads for i in 1:np
            for j in 1:np
                # Compute density of move from s_{i, t} to s_{j, t + 1} * w
                w_i_to_j[i, j] = exp(logpdf_move(xfwd[i, t], xfwd[j, t + 1], zdim, model_move, t, vmap, n_sim, cache) + lw[i])
            end
        end

        # Compute k-sum for each particle j at time t: ∑ₖ f(sⱼ,ₜ₊₁ | sₖ,ₜ) * w
        # * I.e., the summed density of moves from k (t) into each j (t + 1) * w, summed over k
        # * This is a sum over rows (by column) for the w_i_to_j matrix 
        # * It is indexed by j
        ksums = sum(w_i_to_j, dims = 1)

        # Validate ksums
        # * It is possible that the subset of recorded particles at t is not connected to the subset at t + 1
        # * But ksums[j] should not = 0 (below ksums is used for division)
        # * This check guards against ksums = 0 errors 
        if any(ksums .== 0)
            error("The summed density of moves from particle k (t) into j (t + 1) * w, ksums[j], is 0 at time step $t. Boost `n_record` in the filter or the number of particles used in the smoother.")
        end

        # Compute smoothed weights (ws) for each i (t) as a sum over j (t + 1)
        # 1. For each particle j (matrix column) at t + 1
        # 1.1. Divide by jth ksum element 
        # 1.2. Multiply by jth smoothed weight ws element
        # 2. Then sum over columns (j), by row
        # * Note that ksums & ws must be one-row matrix for correct broadcasting 
        # * We return a _vector_ of weights for resampling
        ws = reshape(ws, 1, :) 
        ws = vec(sum((w_i_to_j ./ ksums[:, 1:np]) .* ws[:, 1:np], dims = 2))

        # Validate weights
        if all(x -> x == 0 || isnan(x), ws)
            error("All weights are zero or NaN at time step $t.")
        end
        if !isapprox(sum(ws), 1.0)
            error("Weights do not sum to one at time step $t (sum: $(sum(ws))).")
        end

        # Compute ESS
        ess[t] = 1 / sum(abs2, ws)

        # Resample particles and store for time t using ws
        idx = resample(ws, np)
        xout[:, t] .=  xfwd[idx, t]

        # Update & renormalise weights
        # * ws is the smoothed weight of every particle (i) at time t
        # * This becomes the smoothed every particle (j) at time t + 1
        # * ... before it is updated for every particle (i) at the next time step
        # * Renormalise 
        ws = ws[idx]
        ws = ws ./ sum(ws)

    end

    #### Clean up
    # Reset map_value
    xout = edit_map_value(xout, model_move.map)

    #### Return outputs
    # Follow particle_filter() format
    (
        timesteps   = collect(1:nt), 
        timestamps  = timeline, 
        state       = xout, 
        direction   = nothing, 
        ess         = ess, 
        maxlp       = NaN, 
        convergence = true, 
        trials      = NaN
    )

end