using Distributions
using LogExpFunctions: logistic, log1pexp

export ModelObs
export ModelObsAcousticLogisTrunc
export ModelObsDepthUniformSeabed, ModelObsDepthNormalTruncSeabed
export ModelObsContainer

"""

    Observation models

# `ModelObs`

`ModelObs` is an Abstract Type that groups observation model structures. See below for built-in sub-types. 

# Built-in sub-types 

## `ModelObsAcousticLogisTrunc`

`ModelObsAcousticLogisTrunc` is a `ModelObs` structure for an acoustic observation (0, 1) and a truncated logistic detection probability model. This contains the following fields:
- `sensor_id`: An integer that defines the sensor (receiver) ID;
- `receiver_x`, `receiver_y`: Floats that define the x and y coordinates of the receiver;
- `receiver_alpha`, `receiver_beta`, `receiver_gamma`: Floats that define the parameters of a truncated logistic detection probability model. 

An acoustic observation (``y^{(A)}_{t, k} \\in {0, 1}``) at receiver ``k`` (location ``\\textit{\\textbf{r}}_k = (\\text{receiver\\_x}, \\text{receiver\\_y})``) at time ``t`` is modelled using a Bernoulli probability mass function:

```math
f(y^{(A)}_{t, k} | \\textit{\\textbf{s}}_t) = \\text{Bernoulli}(p_{k,t}(\\textit{\\textbf{s}}_t))
```

where ``p_{k,t}(\\textit{\\textbf{s}}_t)`` is the probability of a detection at receiver ``k`` at time ``t`` given a transmission from location ``\\textit{\\textbf{s}}_t = (x, y)``. `ModelObsAcousticLogisTrunc` represents ``p_{k,t}(\\textit{\\textbf{s}}_t)`` as a logistic function of the Euclidean location between the receiver and the transmitter, according to the equation:
 
```math
p_{k,t}(\\textit{\\textbf{s}}_t) = \\left\\{
\\begin{array}{ll}
(1 + e^{-(\\text{receiver\\_alpha} - \\text{receiver\\_beta} \\cdot |\\textit{\\textbf{s}}_t - \\textit{\\textbf{r}}_k|)})^{-1} & \\text{if } |\\textit{\\textbf{s}}_t - \\textit{\\textbf{r}}_k| < \\text{receiver\\_gamma} \\\\
0 & \\text{otherwise}
\\end{array}
\\right.
```

where ``\\text{receiver\\_gamma}`` is the detection range. 

To simulate an acoustic observation (``y^{(A)}_{t, k} \\in {0, 1}``) from this model, we can draw a sample from a Bernoulli distribution:

```math
y^{(A)}_{t, k} | \\textit{\\textbf{s}}_t \\sim \\text{Bernoulli}(p_{k,t}(\\textit{\\textbf{s}}_t))
```
via `Patter.simulate_obs()`.

## `ModelObsDepthUniformSeabed`

`ModelObsDepthUniformSeabed` is `ModelObs` structure for a depth observation and a uniform depth model. This contains the following fields:

- `sensor_id`: An integer that defines the sensor (tag) ID;
- `depth_shallow_eps`: A float that defines the shallow depth error;
- `depth_deep_eps`: A float that defines the deep depth error;

This model assumes that an individual must be located in an envelope around the seabed, defined by two error terms (`depth_shallow_eps` and `depth_deep_eps`), according to the equation:

```math
f\\left( y_t^{(D)} |  \\textit{\\textbf{s}}_t \\right) =
\\begin{cases} 
z_t & \\text{if } b(\\textit{\\textbf{s}}_t) - \\text{depth\\_shallow\\_eps} \\leq y_t^{(D)} \\leq b(\\textit{\\textbf{s}}_t) + \\text{depth\\_deep\\_eps} \\\\
0 & \\text{otherwise}
\\end{cases}
```

where ``y_t^{(D)}`` is the observed depth, ``b(\\textit{\\textbf{s}}_t)`` is the bathymetric depth in location ``\\textit{\\textbf{s}}_t`` and ``z_t`` is a constant. If `depth_shallow_eps` and `depth_deep_eps` are zero, the individual's depth is required to match the bathymetric depth.

We can simulate observations from this model as follows:

```math
y_t^{(D)} |  \\textit{\\textbf{s}}_t \\sim \\text{Uniform}(b(\\textit{\\textbf{s}}_t) + \\text{depth\\_deep\\_eps}, \\text{min}(b(\\textit{\\textbf{s}}_t) - \\text{depth\\_shallow\\_eps}, 0))
```

via `Patter.simulate_obs()`. If `depth_shallow_eps` and `depth_deep_eps` are set to zero, the `Patter.simulate_obs()` method simply returns the bathymetric depth (`state.map_value`).

## `ModelObsDepthNormalTruncSeabed`

`ModelObsDepthNormalTruncSeabed` is a `ModelObs` structure for a depth observation and a truncated normal model. This contains the following fields:

- `sensor_id`: An integer that defines the sensor (tag) ID;
- `depth_sigma`: A float that defines the standard deviation of the normal distribution;
- `depth_deep_eps`: A float that defines the deep truncation parameter;

This model assumes that an individual must be located in an envelope around the seabed, defined by a normal distribution centred at this location with standard deviation `depth_sigma`: 

```math
f(y_t^{(D)} | \\textit{\\textbf{s}}_t) = \\text{TruncatedNormal}(b(\\textit{\\textbf{s}}_t), \\text{depth\\_sigma}^2, 0, b(\\textit{\\textbf{s}}_t)).
```

We can simulate observations from this model as for previous models via `Patter.simulate_obs()`.

## `ModelObsContainer`

`ModelObsContainer` is a `ModelObs` structure for a container. Containers are a computational device used to mitigate particle degeneracy in the particle filter. Containers define the maximum possible distance of an individual from a location (e.g., receiver) that recorded a future observation. Accordingly, the `ModelObsContainer` structure contains the following fields:
- `sensor_id`: An `integer` that defines the sensor ID (e.g., the receiver ID for the receiver that recorded the next detection);
- `centroid_x`, `centroid_y`: `Float`s that define the x and y coordinates of the container's centroid;
- `radius`: A `Float` that defines the radius of the container; 

An acoustic container defines the region within which an individual must be located according to the receiver(s) at which it was next detected. In the particle filter, particles are permitted or killed depending on whether or not they are compatible with this constraint.

We formalise the acoustic container constraint as follows. Let ``t`` denote the movement timeline (where ``t \\in \\{1,2, \\dots ,T\\}``), ``\\boldsymbol{s}`` denote the individual's latent state (typically (x,y) location), ``\\boldsymbol{y}`` denote acoustic observations (detections, non-detections) at receivers and ``\\boldsymbol{\\theta}`` denote static parameters (in the movement and observation models). The likelihood of the observations ``\\boldsymbol{y}_{1:T}`` given the parameters ``\\boldsymbol{\\theta}`` is expressed as

```math
f(\\boldsymbol{y}_{1:T}\\mid\\boldsymbol{\\theta}) =
\\int
\\prod_{t=1}^{T}
f(\\boldsymbol{y}_{t}\\mid\\boldsymbol{s}_{t},\\boldsymbol{\\theta})
\\, f(\\boldsymbol{s}_{t}\\mid\\boldsymbol{s}_{t-1},\\boldsymbol{\\theta})
\\, d\\boldsymbol{s}_{1:T}.
```

Acoustic containerisation represents a redefinition of the likelihood ``f(\\boldsymbol{y}_{t}\\mid\\boldsymbol{s}_{t},\\boldsymbol{\\theta})``, in a manner that does not affect the likelihood ``f(\\boldsymbol{y}_{1:T}\\mid\\boldsymbol{\\theta})``, with the addition of a characteristic function ``\\chi(t,\\boldsymbol{s}_{t},\\boldsymbol{y}) \\in \\{0,1\\}`` that returns zero or one depending on whether or not the latent state ``\\boldsymbol{s}_{t}`` is compatible with the observations at some future time:

```math
f(\\boldsymbol{y}_{1:T}\\mid\\boldsymbol{\\theta}) =
\\int
\\prod_{t=1}^{T}
f(\\boldsymbol{y}_{t}\\mid\\boldsymbol{s}_{t},\\boldsymbol{\\theta})
\\, \\chi(t,\\boldsymbol{s}_{t},\\boldsymbol{y})
\\, f(\\boldsymbol{s}_{t}\\mid\\boldsymbol{s}_{t-1},\\boldsymbol{\\theta})
\\, d\\boldsymbol{s}_{1:T}.
```

The characteristic function ``\\chi(t,\\boldsymbol{s}_{t},\\boldsymbol{y})`` is defined as follows.  
Let ``\\Delta(t)`` be a function that identifies the next time step (at some time after ``t``) at which detection(s) were recorded.  
At time ``t`` the individual must be within a maximum distance (denoted ``\\text{radius}_{k}``) of each receiver ``k`` that recorded the next detection(s):

```math
\\chi(t,\\boldsymbol{s}_{t},\\boldsymbol{y}) =
\\prod_{k}
\\chi_{k}\\!\\bigl(
\\lVert \\boldsymbol{r}_{k,\\Delta(t)} - \\boldsymbol{s}_{t} \\rVert
\\le \\text{radius}_{k}(\\text{mobility})
\\bigr),
```

where ``\\boldsymbol{r}_{k,\\Delta(t)}`` denotes the location of receiver ``k``.

The radius for receiver ``k`` at time ``t`` depends on the number of time steps from ``t`` until the detection event ``(\\Delta(t)-t)``, a parameter ``\\text{mobility} \\in \\boldsymbol{\\theta}`` that defines the maximum possible moveable distance in between two time steps (i.e., from ``t \\to t+1``) and the receiver’s maximum detection range ``\\gamma \\in \\boldsymbol{\\theta}``:

```math
\\text{radius}(\\text{mobility}) =
(\\Delta(t)-t)\\,\\text{mobility} + \\gamma.
```

In practice, acoustic containers are implemented via the `ModelObsContainer` structure (that is, we treat the container like an additional observation, even though it is simply a redundant use of the data).  `ModelObsContainer` instances must be assembled before the filter run, following the equations above. For relevant time steps, it is necessary to build a `ModelObsContainer` for each of the receiver(s) that recorded a detection at the next detection time step. Each instance contains `centroid_x` and `centroid_y` elements that define the centroid of the container and a `radius` element that defines the maximum possible distance of the individual from the container centroid (see above).

In the particle filter, at each time ``t``, we iterate over all observation‑model structures and update particle weights using a `Patter.logpdf_obs` method. For `ModelObsContainer` instances, the `Patter.logpdf_obs` method computes the distance between the particle ``(i)`` and receiver ``k``’s location (centroid) and adds ``0.0`` or ``\\text{-Inf}`` to the log weight (``\\text{lw}``) depending on whether or not the distance is less than the radius:

```math
\\text{lw}_i =
\\text{lw}_i +
\\log\\!\\bigl(
\\chi_{k}\\!\\bigl(
\\lVert \\boldsymbol{r}_{k,\\Delta(t)} - \\boldsymbol{s}_{t} \\rVert
\\le \\text{radius}_{k}(\\text{mobility})
\\bigr)
\\bigr).
```

This approach is valid because we only kill particles that are incompatible with future observations. A filtering–smoothing algorithm that implements acoustic containerisation thus still approximates ``f(\\boldsymbol{s}_{t}\\mid\\boldsymbol{y}_{1:T})`` once all data have been taken into account.

# Custom sub-types

To define a custom sub-type, such as `ModelObsDepthNormal`, simply define a `struct` that is a sub-type of `Patter.ModelObs`:

```
struct ModelObsDepthNormal <: Patter.ModelObs
    sensor_id::Int64
    depth_sigma::Float64
end
```

For communication with `R`, all sub-types should include a `sensor_id` field. 

Add corresponding methods to simulate observations via `Patter.simulate_obs()` and to evaluate log probabilities via `Patter.logpdf_obs()`. 

# Simulation

`Patter.simulate_obs()` is an internal generic function that simulates observations, given the animal's [`State`](@ref) and a `ModelObs` instance. This accepts the following arguments:
-   `state`: A [`State`](@ref) instance;
-   `model_obs`: A [`ModelObs`](@ref) instance;
-   `t`: An integer that defines the time step;

Methods are implemented for all built-in sub-types. Methods can be defined for new sub-types, such as `ModelObsDepthNormal`, as follows:
```
function Patter.simulate_obs(state::StateCXYZ, model_obs::ModelObsDepthNormal, t::Int64)
    dbn   = truncated(Normal(state.z, model_obs.depth_sigma), 0, state.map_value)
    rand(dbn)
end
```

`Patter.simulate_obs()` is wrapped by [`simulate_yobs()`](@ref) for the simulation of observations.

# Log probabilities 

`Patter.logpdf_obs()` is a generic function that calculates the log probability (density) of an observation, given the animal's [`State`](@ref) and a `ModelObs` instance. This accepts the following arguments:
-   `state`: A `State` instance;
-   `model_obs`: A [`ModelObs`](@ref) instance;
-   `t`: An integer that defines the time step;
-   `obs`: The observation;

Methods are implemented for all built-in sub-types. Methods can be defined for new sub-types, such as `ModelObsDepthNormal`, as follows:
```
function Patter.logpdf_obs(state::State, model_obs::ModelObsDepthNormal, t::Int64, obs::Float64)
    dbn   = truncated(Normal(state.map_value, model_obs.depth_sigma),
                      0.0, state.map_value)
    logpdf(dbn, obs)
  end
```

`Patter.logpdf_obs()` is used in [`particle_filter()`](@ref) to evaluate the log-probability of the data given particle samples.
        
"""
abstract type ModelObs end


#########################
#########################
#### Truncated logistic acoustic observation model

struct ModelObsAcousticLogisTrunc <: ModelObs
    sensor_id::Int64
    # Receiver coordinates
    receiver_x::Float64
    receiver_y::Float64
    # Detection probability parameters
    receiver_alpha::Float64
    receiver_beta::Float64    
    receiver_gamma::Float64
end

function simulate_obs(state::State, model_obs::ModelObsAcousticLogisTrunc, t::Int64)
    # Evaluate the distance between the particle and receiver
    dist = distance(state.x, state.y, model_obs.receiver_x, model_obs.receiver_y)
    # Define probability of detection
    prob = ifelse(dist > model_obs.receiver_gamma, 0.0, logistic(model_obs.receiver_alpha + model_obs.receiver_beta * dist))
    # Define distribution
    rand(Bernoulli(prob)) + 0
end

function logpdf_obs(state::State, model_obs::ModelObsAcousticLogisTrunc, t::Int64, obs::Int64)

    # Evaluate the distance between the particle and receiver
    dist = distance(state.x, state.y, model_obs.receiver_x, model_obs.receiver_y)

    # Calculate log probability given detection (1)
    if obs == 1
        if dist > model_obs.receiver_gamma 
            return -Inf
        else 
            return -log1pexp(-(model_obs.receiver_alpha + model_obs.receiver_beta * dist))
        end 
    
    # Calculate log probability given non detection (0)
    elseif obs == 0
        if dist > model_obs.receiver_gamma
            return 0.0
        else 
            return -log1pexp(model_obs.receiver_alpha + model_obs.receiver_beta * dist)
        end

    else 
      error("Acoustic observations should be coded as 0 or 1.")  
    end

end


#########################
#########################
#### Uniform depth model 

struct ModelObsDepthUniformSeabed <: ModelObs
    sensor_id::Int64
    depth_shallow_eps::Float64
    depth_deep_eps::Float64
end

function simulate_obs(state::State, model_obs::ModelObsDepthUniformSeabed, t::Int64)
    if model_obs.depth_shallow_eps == model_obs.depth_deep_eps == 0.0
        return state.map_value
    end 
    a     = max(0, state.map_value - model_obs.depth_shallow_eps)
    b     = state.map_value + model_obs.depth_deep_eps
    dbn   = Uniform(a, b)
    rand(dbn)
end 

function logpdf_obs(state::State, model_obs::ModelObsDepthUniformSeabed, t::Int64, obs::Float64)
    if model_obs.depth_shallow_eps == model_obs.depth_deep_eps == 0.0
        if abs(state.map_value - obs) < 1.0e-7
            return 0.0
        else
            return -Inf
        end 
    end 
    a     = max(0, state.map_value - model_obs.depth_shallow_eps)
    b     = state.map_value + model_obs.depth_deep_eps
    dbn   = Uniform(a, b)
    logpdf(dbn, obs)
end 


#########################
#########################
#### Truncated normal depth model

struct ModelObsDepthNormalTruncSeabed <: ModelObs 
    sensor_id::Int64
    # Normal distribution variance
    depth_sigma::Float64
    # The individual may be `depth_deep_eps` deeper than the bathymetric depth
    depth_deep_eps::Float64
end 

# Truncated normal depth model log probability 
# * The probability is highest if the individual is on the seabed
# * The individual can be up to model_obs.depth_deep_eps deeper than the seabed
# * Probability decays away from the seabed toward the surface
# * Note the need to normalise logpdfs, since the information content depends on the bathymetric depth
function simulate_obs(state::State, model_obs::ModelObsDepthNormalTruncSeabed, t::Int64)
    dbn   = truncated(Normal(state.map_value, model_obs.depth_sigma), 
                      0.0, state.map_value + model_obs.depth_deep_eps)
    rand(dbn)
end 

function logpdf_obs(state::State, model_obs::ModelObsDepthNormalTruncSeabed, t::Int64, obs::Float64)
    dbn   = truncated(Normal(state.map_value, model_obs.depth_sigma), 
                      0.0, state.map_value + model_obs.depth_deep_eps)
    logpdf(dbn, obs)
end 

#########################
#########################
#### Containers

struct ModelObsContainer <: ModelObs
    sensor_id::Int64
    centroid_x::Float64
    centroid_y::Float64
    radius::Float64
end

# A simulate_obs() method is not currently provided for ModelObsContainer
# * Containers are calculated post-hoc
# function simulate_obs(state::State, model_obs::ModelObsContainer, t::Int64)
# 
# end

function logpdf_obs(state::State, model::ModelObsContainer, t::Int64, obs::Int64)
    # Calculate distance between particle (state) and receiver
    dist = distance(state.x, state.y, model.centroid_x, model.centroid_y)
    # Only particles within model.radius are permitted
    # * radius is a pre-calculated field in model
    # * (For acoustics: radius = receiver_gamma + (receiver_timestep - t) * mobility)
    # return ifelse(dist <= model.radius, log(1.0 / (π * model.radius^2)), -Inf)
    # Return 0.0 or -Inf (CA)
    return ifelse(dist <= model.radius, 0.0, -Inf)
end
