using Base.Threads: @threads
using JLD2
using ProgressMeter: @showprogress

export particle_smoother_two_filter


# batch_indices() helper
# * This determines the UnitRange over which we iterate in _particle_smoother_two_filter() in line with the batch
# * `b` is the index of the current batch
# * `nb` is the total number of batches
# * `nt` is the total number of columns 
function batch_indices(b::Int, nb::Int, nt::Int)
    if nb == 1
        indices = 2:(nt - 1)
    else
        if b == 1
            indices = 2:nt
        else
            if b < nb
                indices = 1:nt
            else
                indices = 1:(nt - 1)
            end
        end
    end
    return indices
end 


# smooth_weights() helper to compute weights
# * Compartmentalising this function massively boosts computation time 
# * `xbwd` and `xfwd` are Vectors of particles from `xbwd[:, t]` and `xfwd[:, t - 1]`
# * `t` is the _actual_ time step 
# * `np` is the number of particles, which is predefined for speed only
# * `t` is the _actual_ time step 
# * other arguments are as defined in particle_smoother_two_filter()
function smooth_weights(xbwd::Vector{<:State}, xfwd::Vector{<:State}, 
                        np::Int,
                        zdim::Bool,
                        model_move::ModelMove,
                        t::Int,
                        vmap::Union{GeoArray, Nothing},
                        n_sim::Int,
                        cache_norm_constants::Union{Dict{<:State,Float64},Nothing})
    w = zeros(np)
    @threads for k in 1:np
        for j in 1:np
            # Evaluate probability density of movement between locations (i.e., the weight) for xbwd[k, t] to xfwd[j, t - 1]
            w[k] += exp(logpdf_move(xbwd[k], xfwd[j], zdim, model_move, t, vmap, n_sim, cache_norm_constants))
        end
    end
    return w
end


# smooth_resample() helper to resample smoothed particles and record ESS
# * `xfwd` and `xbwd` are matrices
# * `t` is the _column index_
# * `np`, `n_fwd_half` and `n_fwd_bwd` are constants (predefined for speed)
function smooth_resample(xfwd::Matrix{<:State}, xbwd::Matrix{<:State}, 
                         w::Vector{Float64}, 
                         t::Int,
                         np::Int, n_fwd_half::Int, n_bwd_half::Int)
    if any(w .> 0)
        # (A) If there are positive weights, normalise, compute ESS & resample
        w .= w ./ sum(w)
        ess_t = 1 / sum(abs2, w)
        # Resample particles from xbwd
        idx    = resample(w, np)
        xout_t = xbwd[idx, t]
    else
        # (B) If all weights are zero:
        # * Sample 50 % of particles from forward filter & 50 % from the backward filter
        # * (Particles are equally weighted thanks to resampling)
        # * Leave ess[t] = NaN (throw warning in outer function)
        ess_t  = NaN
        xout_t = vcat(xfwd[1:n_fwd_half, t], xbwd[1:n_bwd_half, t])
    end
    return xout_t, ess_t
end 


# Internal two-filter smoothing function, wrapped by particle_smoother_two_filter()
# * This function implements the core logic and permits batching in the wrapper function 
# * `timesteps` is a Vector of the time steps for a specific batch e.g., (1, 2, 3 or 4, 5, 6 or 7, 8, 9 10)
# * `indices` is a UnitRange of indices over which we iterate
# - If a single batch, indices is 2:(nt-1) (particles at t = 1 are given by xbwd[:, 1])
# - If batch one of several, indices is 2:nt
# - If batch two of several, indices = 1:nt
# - If batch N of N, indices = 1:(nt - 1) (particles at t = T are given by xfwd[1:, T])
# * `xfwd_init` are the locations from the forward filter
# - This is needed for multi-batch implementations to compute density from xbwd[:, t] to xfwd[:, t-1] at t = 1 if we are on batch ID > 1
# * `xfwd` and `xbwd` are the state (location) matrices (for a selected batch)
# * Other arguments are defined as for the wrapper function

function _particle_smoother_two_filter(; timesteps::Vector{Int}, 
                                         indices::UnitRange{Int},
                                         xfwd_init::Union{Nothing, Vector{<:State}}, 
                                         xfwd::Matrix{<:State}, xbwd::Matrix{<:State}, 
                                         model_move::ModelMove, 
                                         vmap::Union{GeoArray, Nothing} = nothing,
                                         n_sim::Int = 100,
                                         cache::Bool = true)

    #### Check inputs
    size(xfwd) == size(xbwd) || error("Forward and backward sample do not match!")

    #### Set up
    # Identify the dimension of the movement model
    zdim = hasfield(typeof(xfwd[1]), :z)
    # Define smoothed particles matrix 
    # (rows: particles; columns: time steps)
    xout = similar(xfwd)
    # Define particle indices
    # * If the two filters are incompatible (all weights zero), 
    # * ... we sample n_fwd_half and n_bwd_half particles from  
    # * ... the forward and backward filters respectively
    np, nt = size(xout)
    half = (np ÷ 2)
    n_fwd_half = length(1:half)
    n_bwd_half = length((half + 1):np)
    # Initialise ESS vector
    ess = fill(NaN, nt)

    #### Precomputations
    # Precompute normalisation constants if cache = true & n_sim > 0
    # * This is possible for movement models for which the density only depends on `xbwd` fields 
    if cache & n_sim > 0
        cache_norm_constants = logpdf_move_normalisations(xbwd, model_move, vmap, n_sim)
    else
        cache_norm_constants = nothing
    end

    #### Run smoothing for t = 2:(nt - 1), 1:nt or 2:(nt - 1):
    @showprogress desc = "Running two-filter smoother..." for t in indices
        
        # println(t)

        # Compute weights
        if t == 1
            # If the column index t = 1 and we are on batch 1, this step is not activated 
            # (indices = 2:(nt - 1) as we use particles from the backward filter at time = 1)
            # For other batches, we compute the weight from the first column index to the previous position on the forward filter 
            # (This is recorded in xfwd_init)
            w = smooth_weights(xbwd[:, 1], xfwd_init, np, zdim, model_move, timesteps[t], vmap, n_sim, cache_norm_constants)
        else
            w = smooth_weights(xbwd[:, t], xfwd[:, t - 1], np, zdim, model_move, timesteps[t], vmap, n_sim, cache_norm_constants)
        end 

        # Rsample particles & compute ESS using smoothed weights 
        xout[:, t], ess[t] = smooth_resample(xfwd, xbwd, w, t, np, n_fwd_half, n_bwd_half)

    end 

    # Return xout and ess
    return (xout = xout, ess = ess)

end

"""
    particle_smoother_two_filter(; timeline::Vector{DateTime}, 
                                   xfwd::Union{Matrix{<:State}, Vector{String}}, 
                                   xbwd::Union{Matrix{<:State}, Vector{String}}, 
                                   model_move::ModelMove, 
                                   vmap::Union{GeoArray, Nothing} = nothing, 
                                   n_sim::Int = 100, 
                                   cache::Bool = true, 
                                   batch::Union{Nothing, Vector{String}} = nothing)

A two-filter particle smoother that samples from `f(s_t | y_{1:T})` for `t ∈ 1:T`.

# Arguments (keywords)

- `timeline`: A `Vector{DateTime}` of ordered, regularly spaced time stamps that defines the time steps for the simulation;
- `xfwd`, `xbwd`: Particles from the forward and backward filters (see [`particle_filter()`](@ref)), supplied as:
    - A `Matrix` of [`State`](@ref)s (in memory);
    - A `Vector` of file paths, if [`particle_filter()`](@ref) was implemented with `batch`;
- `model_move`: A [`ModelMove`](@ref) instance;
- `vmap`: (optional) A `GeoArray` that defines the 'validity map' (see [`Patter.logpdf_move()`](@ref));
- `n_sim`: An integer that defines the number of Monte Carlo simulations (see [`Patter.logpdf_move()`](@ref));
- `cache`: A `Bool` that defines whether or not to precompute and cache movement density normalisation constants (see [`Patter.logpdf_move()`](@ref));
- (optional) `batch`: A `Vector` of `.jld2` file paths for particles (see [`particle_filter()`](@ref));

# Details

[`particle_smoother_two_filter()`](@ref) smooths particles from the particle filter (see [`particle_filter()`](@ref)). The `timeline` from the particle filter should be supplied as well as a `Matrix` of particles from a forward run and a backward run (or a `Vector` of file paths to those matrices). The two filter smoother works by iteratively resampling particles in line with the probability density of movement between particles from the backward filter at time `t` and particles from the forward filter at time `t - 1`. [`Patter.logpdf_move()`](@ref) is an internal function that evaluates the log probability of a movement step between particles. This function wraps the [`Patter.logpdf_step()`](@ref) generic. Methods are provided for built-in [`State`](@ref) and [`ModelMove`](@ref) sub-types. To use custom sub-types, a corresponding [`Patter.logpdf_step()`](@ref) method should be provided. In [`particle_smoother_two_filter()`](@ref), the `vmap` and `n_sim` arguments support the calculate of probability densities (see [`Patter.logpdf_move()`](@ref)). For movement models for which the density only depends on fields in `xbwd` and `xfwd`, set `cache = true` to precompute and store normalisation constants for density calculations for unique `xbwd` elements. Note that since typically only a subsample of particles from [`particle_filter()`](@ref) are retained in memory, it is not guaranteed that valid moves will exist between particle pairs at all time steps. At time step(s) in which the two filters are incompatible, 50 % of particles are retained from the forward filter and 50 % from the backward filter with a warning. The effective sample size at such time steps is set to `NaN`, providing an index and counter for problematic time steps (see Returns). Batching is only implemented if the inputs (`xfwd`, and `xbwd`) and outputs (via `batch`) are batched (and contain the same number of batches). 

# Returns 

- A [`Particles`](@ref) structure;

# See also

- [`particle_filter()`](@ref) implements the particle filter;
- [`Patter.logpdf_step()`](@ref), [`logpdf_move_normalisation()`](@ref) and [`Patter.logpdf_move()`](@ref) evaluate the log probability (density) of movement between two [`State`](@ref)s;
- [`particle_smoother_two_filter()`](@ref) implements the two-filter particle smoother;

# Source 

Fearnhead, P., Wyncoll, D., Tawn, J., [2010](https://doi.org/10.1093/biomet/asq013). A sequential smoothing algorithm with linear computational cost. Biometrika 97, 447–464.

"""
function particle_smoother_two_filter(; timeline::Vector{DateTime}, 
                                        xfwd::Union{Matrix{<:State}, Vector{String}}, xbwd::Union{Matrix{<:State}, Vector{String}}, model_move::ModelMove, 
                                        vmap::Union{GeoArray, Nothing} = nothing, 
                                        n_particle::Union{Nothing, Int} = nothing,
                                        n_sim::Int = 100, 
                                        cache::Bool = true, 
                                        batch::Union{Nothing, Vector{String}} = nothing)

    #### Initialise
    call_start = now()

    #### Prepare batch settings
    # Check xfwd & xbwd sizes (of matrices or files) match
    size(xfwd) == size(xbwd) || error("Forward and backward sample do not match!")
    # We start by assuming 1 batch and set do_batch to false
    nb = 1
    do_batch = false
    # Handle batching
    if (isa(xfwd, Vector{String}))
        # If inputs are batched, we enforce the same number of output batches
        !isnothing(batch) || error("Batched inputs required batched outputs.")
        @assert length(xfwd) == length(xbwd) == length(batch)
        # Record batch file paths
        xfwd_batch = deepcopy(xfwd)
        xbwd_batch = deepcopy(xbwd)
        # Update number of batches and do_batch 
        nb         = length(xfwd_batch)
        do_batch   = true
    else
        isnothing(batch) || error("`batch` is only implemented if inputs are batched.")
    end 

    #### Prepare smoother objects
    # Initialise state matrices/vectors
    xout      = nothing
    xfwd_init = nothing
    # Initialise  ESS vector (over full timetime)
    nt   = length(timeline)
    np   = NaN              # read np below 
    ess  = fill(NaN, nt)    # zeros(nt)
    # Define indices for batches
    # * timesteps_by_batch is a Vector of time steps (one element per batch), as defined in the filter
    # - This is used to update the global ess Vector
    # - E.g., 1, 2, 3,    4, 5, 6,    7, 8, 9, 10
    timesteps_by_batch = split_indices(collect(1:length(timeline)), nb)

    #### Run smoothing (over batches as required)
    for b in 1:nb

        # Read batch if required
        if do_batch
            @load xfwd_batch[b] xfwd
            @load xbwd_batch[b] xbwd
        end 

        # (optional) Select particles 
        if !isnothing(n_particle)
            xfwd = xfwd[1:n_particle, :]
            xbwd = xbwd[1:n_particle, :]
        end 

        println(xfwd_batch[b])
        println(xbwd_batch[b])
        println(size(xfwd))
        println(size(xbwd))

        # Define indicies for batch
        timesteps = timesteps_by_batch[b]
        indices   = batch_indices(b, nb, size(xfwd, 2))

        # Run smoother (for batch)
        # If b = 1, xfwd_init = nothing as xfwd_init is not required 
        bout = _particle_smoother_two_filter(timesteps  = timesteps, 
                                             indices    = indices,
                                             xfwd_init  = xfwd_init,
                                             xfwd       = xfwd, 
                                             xbwd       = xbwd,
                                             model_move = model_move,
                                             vmap       = vmap,
                                             n_sim      = n_sim,
                                             cache      = cache)

        # Update particles for t = 1 and t = T 
        if b == 1
            np = size(xbwd, 1)            # Record np once here for convenience (used below)
            bout.xout[:, 1] = xbwd[:, 1]  # Update particles
        end
        if b == nb
            bout.xout[:, end] = xfwd[:, end]
        end

        # Record ess
        # (ESS at t = 1 and t = T is updated outside of the loop over batches)
        ess[timesteps] = bout.ess

        # Record particles 
        if do_batch
            # Record xfwd_init
            # The last locations for this batch are the locations at t - 1 on the next batch
            xfwd_init = xfwd[:, end]
            # Write particles to file and leave xout = nothing 
            @save batch[b] xsmo = bout.xout
        else 
            xout = bout.xout
        end 

    end 

    #### Update ESS
    # Update ESS for t = 1 and t = T 
    ess[1]  = Float64(np)  # = 1 / sum(abs2, w) given equal weights from filter
    ess[nt] = Float64(np)
    # Evaluate smoothing success
    # * Set convergence = false if > 5 % NaN ESS
    nan_count   = count(isnan, ess)
    nan_perc    = round(nan_count / length(ess) * 100, digits=2)
    convergence = nan_perc <= 5 ? true : false
    if nan_count > 0
        julia_warning("All smoothing weights (from xbwd[k, t] to xfwd[j, t - 1]) are zero at $nan_count time step(s) ($nan_perc %).")
    end

    #### Return outputs
    particulate("smoother: two-filter", call_start, 
                timeline, xout, ess, fill(NaN, length(timeline)), 
                np, NaN, convergence)

end