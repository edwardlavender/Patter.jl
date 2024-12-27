using Rasters, LibGEOS
using DataFrames, DataFramesMeta
using Dates

export simulate_states_init

"""
    simulate_states_init(; map::Rasters.Raster,
                           timeline::Vector{DateTime}, 
                           state_type::Type{<: State}, 
                           xinit, 
                           model_move::ModelMove, 
                           datasets, 
                           model_obs_types, 
                           n_particle::Int, 
                           direction::String = "forward", 
                           output = "DataFrame")

Simulate a DataFrame or Vector of initial [`State`](@ref)s for the simulation of movement paths in [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref).

# Arguments

- `map`: A `Rasters.Raster` that defines the study area for the simulation. Here, `map` is used to:
  - Sample initial coordinates, via `Patter.coords_init()`, if `xinit = nothing`;

- `timeline`, `model_move`, `datasets`, `model_obs_types`, `direction`: Additional arguments used to restrict `map`, via `Patter.map_init()`, before sampling initial states.
  - `timeline`: A sorted, `DateTime` vector of regularly spaced time stamps that defines the timeline for the simulation;
  - `model_move`: A [`ModelMove`](@ref) instance;
  - (optional) `datasets`: A `Vector` of observation datasets or `nothing`;
  - (optional) `model_obs_types`: A `Vector` of [`ModelObs`](@ref) sub-types or `nothing`;
  - `direction`: A `String` string that defines the direction of the simulation (`"forward"` or `"backward"`);

- `state_type`: The [`State`](@ref) sub-type. Here, `state` is used to:
  - Convert sampled coordinates to initial states, via `Patter.states_init()`, if `xinit = nothing`;

- `xinit`: (optional) A `DataFrame` of initial states, with one column for each state dimension.

- `n_particle`: An `integer` that defines the number of simulated states:
  - If `xinit = nothing`, `n_particle` specifies the number of simulated states via `Patter.coords_init()`;
  - If `xinit` is supplied but there are not `n_particle` initial states, `n_particle` initial states are re-sampled from `xinit` with replacement;

- `output`: A `String` that defines the output format:
  - `"DataFrame"` returns a `DataFrame` (for communication with `R`);
  - `"Vector"` returns a Vector of [`State](@ref)(s);

# Details

These functions support the simulation of initial states for animal movement walks in [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref).

If `xinit = nothing`, initial coordinates are sampled from `map`.

The region(s) within `map` from which initial coordinates are sampled can be optionally restricted by the provision of the observation datasets and the associated model sub-types (via `Patter.map_init_iter`). This option does not apply to [`simulate_path_walk()`](@ref) but is used in [`particle_filter()`](@ref) where observation models are used. In this instance, `Patter.map_init_iter` iterates over each model and uses the `Patter.map_init()` method to update `map`. The following methods are implemented:
  - Default. The default method returns `map` unchanged.
  - `model_obs_type::ModelObsAcousticLogisTrunc`. This method uses acoustic observations to restrict `map` via Lavender et al.'s ([2023](https://doi.org/10.1111/2041-210X.14193)) acoustic--container algorithm. The function identifies the receiver(s) that recorded detection(s) immediately before, at and following the first time step (`timeline[start]`, where `start` is `1` if `direction = "forward"` and `length(timeline)` otherwise). The 'container' within which the individual must be located from the perspective of each receiver is defined by the time difference and the individual's mobility (that is, the maximum moveable distance the individual could move between two time steps), which must be specified in `model_move.mobility`. The intersection between all containers defines the possible locations of the individual at the first time step.
  - `model_obs_type::ModelObsDepthUniform`. This method uses the depth observations to restrict `map` (which should represent the bathymetry in a region). The individual must be within a region in which the observed depth at `timeline[start]` is within a depth envelope around the bathymetric depth defined by the parameters `depth_shallow_eps` and `depth_deep_eps` (see [`ModelObs`](@ref)). (If there is no observation at `timeline[start]`, `map` is returned unchanged.)
  - `model_obs_type::ModelObsDepthNormalTrunc`. This method also uses depth observations to restrict `map`. The individual must be in a location where the bathymetric depth plus the `depth_deep_eps` parameter at `timeline[start]` is greater than or equal to the observed depth at `timeline[start]` (see [`ModelObs`](@ref)). (If there is no observation at `timeline[start]`, `map` is returned unchanged.)

To handle custom [`ModelObs`](@ref) sub-types, process `map` beforehand or write an appropriate `Patter.map_init()` method.

Using `map`, a `DataFrame` of `n_particle` initial coordinates (`map_value`, `x`, `y`) is sampled using `Patter.coords_init()`. Additional state dimensions are added, as required depending on the `state_type`, via a `Patter.states_init()` method. For custom [`State`](@ref) sub-types, a corresponding `Patter.states_init()` method is required (or supply `xinit` yourself).

If `xinit()` is provided and `n_particle` initial states are provided, `xinit` is returned unchanged. Otherwise, `n_particle` initial states are resampled from `xinit`, with replacement, and returned.

# Returns

[`simulate_states_init()`](@ref) returns:
  - If `output = "DataFrame"`, a `DataFrame`, with `n_particle` rows, and one column for each [`State`](@ref) dimension;
  - If `output = "Vector"`, a `Vector` of [`State`](@ref) instances, with `n_particle` elements;

# See also

These functions are used to initialise simulated movement trajectories in [`simulate_path_walk()`](@ref) and [`particle_filter()`](@ref). See also:

* [`State`](@ref) and [`ModelMove`](@ref) for [`State`](@ref) and movement model sub-types;
* [`Patter.simulate_step()`](@ref) and [`Patter.simulate_move()`](@ref) to simulate new [`State`](@ref)s;
* [`simulate_path_walk()`](@ref) to simulate animal movement paths (via [`ModelMove`](@ref));
* [`simulate_yobs()`](@ref) to simulate observations arising from simulated movements (via [`ModelObs`](@ref));

"""
function simulate_states_init end 

# The default method (for unknown models) returns `map`
function map_init(map::Rasters.Raster, 
                  timeline::Vector{DateTime}, 
                  model_move::ModelMove, 
                  dataset,                        # ::DataFrame,
                  model_obs_type,                 # ::Type{<: ModelObs}, 
                  direction::String = "forward")
    return map
end

# For ModelObsAcousticLogisTrunc, we use the AC algorithm
function map_init(map::Rasters.Raster, 
                  timeline::Vector{DateTime}, 
                  model_move::ModelMove, 
                  dataset::DataFrame,
                  model_obs_type::Type{ModelObsAcousticLogisTrunc}, 
                  direction::String = "forward")

    #### Check user inputs
    dataset = deepcopy(dataset)
    check_names(dataset, "receiver_gamma")

    #### Define a list of container info
    # Define a timeline of detections 
    t1 = ifelse(direction == "forward", first(timeline), last(timeline))
    secs = diffsecs(timeline[2], timeline[1])
    @chain dataset begin
        @subset! :obs .== 1                    
        sort!(:timestamp) 
        @transform! :gap = diffsecs(:timestamp, t1) ./ secs
    end
    if nrow(dataset) == 0
        return map
    end 
    # Define the first detection(s) before t1
    det_before = dataset[dataset.gap .< 0, :]
    if nrow(det_before) > 0 
        det_before = det_before[det_before.gap .== maximum(det_before.gap), :]
    end 
    # Define detection(s) at t1
    det_at = dataset[dataset.gap .== 0, :]
    # Define the first detection(s) after t1
    det_after = dataset[dataset.gap .> 0, :]
    if nrow(det_after) > 0 
        det_after = det_after[det_after.gap .== minimum(det_after.gap), :]
    end
    # Define container info
    cinfo = vcat(det_before, det_at, det_after)
    @transform! cinfo :buffer = :receiver_gamma + (abs.(:gap) * model_move.mobility)

    #### Define acoustic containers 
    # TO DO 
    # * Exclude regions overlapping with receivers that did not record detections
    # Define container as SpatVector
    vcontainers = [LibGEOS.buffer(LibGEOS.Point(row.receiver_x, row.receiver_y),
                   row.buffer, 1000) for row in eachrow(cinfo)]
    vcontainer = spatIntersect(vcontainers)
    # Define container as SpatRaster
    # * Mask map
    # * mask() returns missing (outside container) & original values
    # * mask() also retains missing on land (desirable)
    # * mask!() converts returns true/false (undesirable)
    map = mask(map, with = vcontainer)
    return map

end 

# For ModelObsDepthUniform, we restrict map using the depth observation
function map_init(map::Rasters.Raster, 
    timeline::Vector{DateTime}, 
    model_move::ModelMove, 
    dataset::DataFrame,
    model_obs_type::Type{ModelObsDepthUniform}, 
    direction::String = "forward")

    # Check dataset
    check_names(dataset, ["timestamp", "obs", "depth_shallow_eps", "depth_deep_eps"])
    # Identify the first depth observation 
    t1  = ifelse(direction == "forward", first(timeline), last(timeline))
    pos = findall(dataset.timestamp .== t1)
    if length(pos) == 0
        return map
    elseif length(pos) != 1
        error("Multiple depth observations recorded at the first time stamp.")
    end 
    depth = dataset.obs[pos]
    # Define the corresponding structure parameters
    depth_shallow_eps = dataset.depth_shallow_eps[pos]
    depth_deep_eps    = dataset.depth_deep_eps[pos]
    # Mask map between limits (false/true -> NaN/true)
    # * map + depth_deep_eps must be >= depth
    msk = (map .- depth_shallow_eps .<= depth) .& (map .+ depth_deep_eps .>= depth)
    msk = classify(msk, false => missingval(map))
    map = mask(map, with = msk)
    return map 

end 

# For ModelObsDepthNormalTrunc, we restrict map using the depth observation
function map_init(map::Rasters.Raster, 
    timeline::Vector{DateTime}, 
    model_move::ModelMove, 
    dataset::DataFrame,
    model_obs_type::Type{ModelObsDepthNormalTrunc}, 
    direction::String = "forward")
    # Check dataset
    check_names(dataset, ["timestamp", "obs", "depth_sigma", "depth_deep_eps"])
    # Identify the first depth observation 
    t1  = ifelse(direction == "forward", first(timeline), last(timeline))
    pos = findall(dataset.timestamp .== t1)
    if length(pos) == 0
        return map
    elseif length(pos) != 1
        error("Multiple depth observations recorded at the first time stamp.")
    end 
    depth = dataset.obs[pos]
    # Define the corresponding structure parameters
    depth_deep_eps = dataset.depth_deep_eps[pos]
    # Mask map between limits (false/true -> NaN/true)
    # * map + depth_deep_eps must be >= depth
    msk = (map .+ depth_deep_eps .>= depth)
    msk = classify(msk, false => missingval(map))
    map = mask(map, with = msk)
    return map

end 

# Iteratively update map (according for each input dataset) using map_init methods
function map_init_iter(
    map::Rasters.Raster, 
    timeline::Vector{DateTime}, 
    model_move::ModelMove, 
    datasets,                # ::Vector,
    model_obs_types,         # ::Vector{Type{<: ModelObs}}, 
    direction::String = "forward")

    if isnothing(datasets)
        return(map)
    end 

    for i in eachindex(model_obs_types)
        map = map_init(map, timeline, model_move, datasets[i], model_obs_types[i], direction)
    end 
    return map

end 

# Sample .n initial coordinates (map_value, x, y) from a SpatRaster
function coords_init(map, size) 
    # Sample initial coordinates (x, y, z)
    xinit = spatSample(x = map, size = size)
    if nrow(xinit) != size
        error("sime_coord_init() failed to generate $size samples.")
    end 
    return xinit
end

# Convert a `DataFrame` of coordinates (map_value, x, y) to a DataFrame with all state dimensions
function states_init(state_type::Type{<: State}, coords)
 error("For custom states, you need to define a `Patter.states_init()` method or provide `.xinit`.")
end 

function states_init(state_type::Type{StateXY}, coords)
    return coords
end 

function states_init(state_type::Type{StateXYZ}, coords)
    N = nrow(coords)
    # Add z coordinate
    @transform! coords :z = :map_value .* rand(N)
    return coords
end 

function states_init(state_type::Type{StateCXY}, coords)
    N = nrow(coords)
    # Add heading
    @transform! coords :heading = rand(N) .* 2 .* π
    return coords
end 

function states_init(state_type::Type{StateCXYZ}, coords)
    N = nrow(coords)
    # Add z coordinate
    @transform! coords :z = :map_value .* rand(N)
    # Add heading
    @transform! coords :heading = rand(N) .* 2 .* π
    return coords
end 

# Simulate states
function simulate_states_init(;
    map::Rasters.Raster, 
    timeline::Vector{DateTime},
    state_type::Type{<: State},
    xinit,
    model_move::ModelMove, 
    datasets,                    # ::Vector, 
    model_obs_types,             # ::Vector, 
    n_particle::Int,
    direction::String = "forward", 
    output = "DataFrame")

    if !(output in ["DataFrame", "Vector"])
        error("`output` must be 'DataFrame' or 'Vector'.")
    end

    #### If un-provided, sample `.xinit`
    if isnothing(xinit) 
        
        #### Define an initial map from which to sample
        # We use the observation datasets to restrict (if possible) the input `map` for sampling
        # We do not need to deepcopy map b/c it is never updated by reference
        # But, if this approach changes, a deepcopy is required for safe application
        # of the algorithms from R for multiple time series dependent on the same initial map (`env_init`)
        if !isnothing(model_obs_types) 
            map = map_init_iter(map, 
                                timeline,
                                model_move,
                                datasets,
                                model_obs_types,
                                direction)
        end 

        #### Simulate initial states
        # Define initial coordinates (map_value, x, y)
        coords = coords_init(map, n_particle)
        # Add additional state dimensions
        # * User defined State structures require:
        # - states_init method or `.xinit`
        xinit = states_init(state_type, coords)

    end 
    
    #### Re-sample `xinit` `n_particle` times (if required)
    # This is required if:
    # * The user supplies `.xinit` with fewer than specified rows
    if nrow(xinit) != n_particle
        xinit = xinit[sample(1:nrow(xinit), n_particle, replace = true), :]
    end

    #### Return outputs
    if output == "DataFrame"
        return xinit
    elseif output == "Vector"
        return julia_get_xinit(state_type, xinit)
    end 

end 

