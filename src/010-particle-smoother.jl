using Base.Threads: @threads
using ProgressMeter: @showprogress

export two_filter_smoother

"""
    two_filter_smoother(; timeline::Vector{DateTime}, 
                        xfwd::Matrix, xbwd::Matrix, 
                        model_move::ModelMove, 
                        vmap::Union{GeoArray, Nothing}, 
                        n_sim::Int, 
                        cache::Bool)

A two-filter particle smoother that samples from `f(X_t | {Y_1 ... Y_T}) for t ∈ 1:T`.

# Arguments (keywords)

- `timeline`: A `Vector{DateTime}` of ordered, regularly spaced time stamps that defines the time steps for the simulation;
- `xfwd`: A `Matrix` of [`State`](@ref)s from the forward filter (see [`particle_filter()`](@ref));
- `xbwd`: A `Matrix` of [`State`](@ref)s from the backward filter (see [`particle_filter()`](@ref));
- `model_move`: A [`ModelMove`](@ref) instance;
- `vmap`: (optional) A `GeoArray` that defines the 'validity map' (see [`logpdf_move()`](@ref));
- `n_sim`: An integer that defines the number of Monte Carlo simulations (see [`logpdf_move()`](@ref));
- `cache`: A `Bool` that defines whether or not to precompute and cache movement density normalisation constants (see [`logpdf_move()`](@ref));

# Details

[`two_filter_smoother()`](@ref) smooths particles from the particle filter (see [`particle_filter()`](@ref)). The `timeline` from the particle filter should be supplied as well as a `Matrix` of particles from a forward run and a backward run. The two filter smoother works by iteratively resampling particles in line with the probability density of movement between particles from the backward filter at time `t` and particles from the forward filter at time `t - 1`. [`logpdf_move()`](@ref) is an internal function that evaluates the log probability of a movement step between particles. This function wraps the [`logpdf_step()`](@ref) generic. Methods are provided for built-in [`State`](@ref) and [`ModelMove`](@ref) sub-types. To use custom sub-types, a corresponding [`logpdf_step()`](@ref) method should be provided. In [`two_filter_smoother()`](@ref), the `vmap` and `n_sim` arguments support the calculate of probability densities (see [`logpdf_move()`](@ref)). For movement models for which the density only depends on fields in `xbwd` and `xfwd`, set `cache = true` to precompute and store normalisation constants for density calculations for unique `xbwd` elements. Note that since typically only a subsample of particles from [`particle_filter()`](@ref) are retained in memory, it is not guaranteed that valid moves will exist between particle pairs at all time steps. At time step(s) in which the two filters are incompatible, 50 % of particles are retained from the forward filter and 50 % from the backward filter with a warning. The effective sample size at such time steps is set to `NaN`, providing an index and counter for problematic time steps (see Returns). 

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
function two_filter_smoother(;timeline::Vector{DateTime}, 
                             xfwd::Matrix, xbwd::Matrix, model_move::ModelMove, 
                             vmap::Union{GeoArray, Nothing} = nothing, 
                             n_sim::Int = 100, 
                             cache::Bool = true)

    #### Check inputs
    size(xfwd) == size(xbwd) || error("Forward and backward sample do not match!")

    #### Set up
    # Identify dimension of input state
    # * This is used to check if states are valid (incl. for normalisation simulation)
    zdim = hasfield(typeof(xfwd[1]), :z)
    # Define smoothed particles matrix 
    # (rows: particles; columns: time steps)
    xout = similar(xfwd)
    xout[:, 1] = xbwd[:, 1]
    xout[:, end] = xfwd[:, end]
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

    #### Precomputations
    # Precompute normalisation constants if cache = true
    # * This is possible for movement models for which the density only depends on `xbwd` fields 
    if cache
        cache_norm_constants = logpdf_move_normalisations(xbwd, model_move, vmap, n_sim)
    else
        cache_norm_constants = nothing 
    end  

    #### Run smoothing
    @showprogress desc = "Running two-filter smoother..." for t in 2:(nt - 1)

        # Compute weights
        # * Multi-threaded implementation is somewhat faster with Dict() cache_norm_constants formulation 
        w = zeros(np)
        @threads for k in 1:np
            for j in 1:np
                # Evaluate probability density of movement between locations (i.e., the weight)
                w[k] += exp(logpdf_move(xbwd[k, t], xfwd[j, t - 1], zdim, model_move, t, vmap, n_sim, cache_norm_constants))
             end
        end
        
        # Validate weights & implement resampling 
        if any(w .> 0)
            # (A) If there are positive weights, normalise, compute ESS & resample
            w .= w ./ sum(w)
            ess[t] = 1 / sum(abs2, w)
            # Resample particles from xbwd & store
            idx = resample(w, np)
            xout[:, t] =  xbwd[idx, t]
        else 
            # (B) If all weights are zero, set a warning flag 
            # * For speed, the warning is thrown only once outside the for loop 
            # * (as it may occur at multiple time steps)
            warn_zero_weights = true
            # Set ESS[t] = NaN
            ess[t] = NaN
            # Sample 50 % of particles from forward filter
            # * We simply select the first n_fwd_half particles since
            # * ... particles are equally weighted thanks to resampling  
            xout[1:half , t] = xfwd[1:n_fwd_half, t]
            # Sample 50 % of particles from backward filter 
            xout[(half + 1):end, t] = xbwd[1:n_bwd_half, t] 

        end 

    end

    #### Implement post-processing 
    # Warn for smoothing failures
    if warn_zero_weights
        nan_count = count(isnan, ess)
        nan_perc  = round(nan_count / length(ess) * 100, digits = 2)
        julia_warning("All smoothing weights (from xbwd[k, t] to xfwd[j, t - 1]) are zero at $nan_count time step(s) ($nan_perc %).")
    end 

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