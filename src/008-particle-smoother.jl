using Base.Threads: @threads
using ProgressMeter: @showprogress

export two_filter_smoother

"""
    two_filter_smoother(; timeline::Vector{DateTime}, xfwd::Matrix, xbwd::Matrix, model_move::ModelMove, box, n_sim::Int)

A two-filter particle smoother that samples from `f(X_t | {Y_1 ... Y_T}) for t ∈ 1:T`.

# Arguments (keywords)

- `timeline`: A `Vector{DateTime}` of ordered, regularly spaced time stamps that defines the time steps for the simulation;
- `xfwd`: A `Matrix` of [`State`](@ref)s from the forward filter (see [`particle_filter()`](@ref));
- `xbwd`: A `Matrix` of [`State`](@ref)s from the backward filter (see [`particle_filter()`](@ref));
- `model_move`: A [`ModelMove`](@ref) instance;
- `box`: (optional) A `NamedTuple` (`min_x`, `max_x`, `min_y`, `max_y`) that defines a 'mobility box' (see [`logpdf_move()`](@ref));
- `n_sim`: An integer that defines the number of Monte Carlo simulations (see [`logpdf_move()`](@ref));

# Details

[`two_filter_smoother()`](@ref) smooths particles from the particle filter (see [`particle_filter()`](@ref)). The `timeline` from the particle filter should be supplied as well as a `Matrix` of particles from a forward run and a backward run. The two filter smoother works by iteratively resampling particles in line with the probability density of movement between particles from the backward filter at time `t` and particles from the forward filter at time `t - 1`. [`logpdf_move()`](@ref) is an internal function that evaluates the log probability of a movement step between particles. This function wraps the [`logpdf_step()`](@ref) generic. Methods are provided for built-in [`State`](@ref) and [`ModelMove`](@ref) sub-types. To use custom sub-types, a corresponding [`logpdf_step()`](@ref) method should be provided. In [`two_filter_smoother()`](@ref), the `box` and `n_sim` arguments support the calculate of probability densities (see [`logpdf_move()`](@ref)). 

# Returns 

- A `NamedTuple`, of the same format as returned by [`particle_filter()`](@ref), with the following fields:
    - `timesteps`
    - `timestamps`
    - `direction`: `nothing`
    - `state`
    - `ess`
    - `maxlp`: `NaN`
    - `convergence`: `true`

# See also

- [`particle_filter()`](@ref) implements the particle filter;
- [`logpdf_step()`](@ref), [`logpdf_move_normalisation()`](@ref) and [`logpdf_move()`](@ref) evaluate the log probability (density) of movement between two [`State`](@ref)s;
- [`two_filter_smoother()`](@ref) implements the two-filter particle smoother;

# Source 

Fearnhead, P., Wyncoll, D., Tawn, J., [2010](https://doi.org/10.1093/biomet/asq013). A sequential smoothing algorithm with linear computational cost. Biometrika 97, 447–464.

"""
function two_filter_smoother(;timeline::Vector{DateTime}, xfwd::Matrix, xbwd::Matrix, model_move::ModelMove, box = nothing, n_sim::Int = 100)

    #### Check inputs
    size(xfwd) == size(xbwd) || error("Forward and backward sample do not match!")

    #### Set up
    # Dimension of input state
    zdim = hasfield(typeof(xfwd[1]), :z)
    # Matrix for smoothed particles
    xout = similar(xfwd)
    xout[:, 1] .= xbwd[:, 1]
    xout[:, end] .= xfwd[:, end]
    np, nt = size(xout)
    # ESS vector
    ess = zeros(nt)
    ess[1] = NaN
    ess[nt] = NaN 
    # Least recently used cache 
    cache = LRU{eltype(xfwd), Float64}(maxsize = np)

    #### Run two-filter formula
    @showprogress desc = "Running two-filter smoother..." for t in 2:(nt - 1)
        w = zeros(np)
        @threads for k in 1:np
            for j in 1:np
                # Evaluate probability density of movement between locations (i.e., the weight)
                w[k] += exp(logpdf_move(xbwd[k, t], xfwd[j, t - 1], zdim, model_move, t, box, n_sim, cache))
            end
        end
        # Normalise weights & evaluate ESS
        w .= w ./ sum(w)
        ess[t] = 1 / sum(abs2, w)
        # Resample particles & store
        idx = resample(w, np)
        xout[:, t] .=  xbwd[idx, t]
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
        convergence = true
    )

end