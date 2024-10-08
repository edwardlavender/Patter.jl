using Distributions
using LogExpFunctions: logistic, log1pexp

export ModelObs
export ModelObsAcousticLogisTrunc, ModelObsAcousticLogisStd
export ModelObsDepthUniform, ModelObsDepthNormalTrunc, ModelObsDepthNormal


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
via [`Patter.simulate_obs()`](@ref).

## `ModelObsDepthUniform`

`ModelObsDepthUniform` is `ModelObs` structure for a depth observation and a uniform depth model. This contains the following fields:

- `sensor_id`: An integer that defines the sensor (tag) ID;
- `depth_shallow_eps`: A float that defines the shallow depth error;
- `depth_deep_eps`: A float that defines the deep depth error;

This model assumes that an individual must be located in an envelope around the bathymetric depth, defined by two error terms (`depth_shallow_eps` and `depth_shallow_eps`), according to the equation:

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

via [`Patter.simulate_obs()`](@ref). If `depth_shallow_eps` and `depth_deep_eps` are set to zero, the [`Patter.simulate_obs()`](@ref) method simply returns the bathymetric depth (`state.map_value`).

## `ModelObsDepthNormalTrunc`

`ModelObsDepthNormalTrunc` is a `ModelObs` structure for a depth observation and a truncated normal model. This contains the following fields:

- `sensor_id`: An integer that defines the sensor (tag) ID;
- `depth_sigma`: A float that defines the standard deviation of the normal distribution;
- `depth_deep_eps`: A float that defines the deep truncation parameter;

This model assumes that an individual must be located in an envelope around the bathymetric depth, defined by a normal distribution centred at this location with standard deviation `depth_sigma`: 

```math
f(y_t^{(D)} | \\textit{\\textbf{s}}_t) = \\text{TruncatedNormal}(b(\\textit{\\textbf{s}}_t), \\text{depth\\_sigma}^2, 0, b(\\textit{\\textbf{s}}_t)).
```

We can simulate observations from this model as for previous models via [`Patter.simulate_obs()`](@ref).

# Custom sub-types

To define a custom sub-type, such as `ModelObsDepthNormal`, simply define a `struct` that is a sub-type of `Patter.ModelObs`:

```
struct ModelObsDepthNormal <: Patter.ModelObs
    sensor_id::Int64
    depth_sigma::Float64
end
```

For communication with `R`, all sub-types should include a `sensor_id` field. 

Add corresponding methods to simulate observations via [`Patter.simulate_obs()`](@ref) and to evaluate log probabilities via [`Patter.logpdf_obs()`](@ref). 

# Simulation

[`Patter.simulate_obs()`](@ref) is an internal generic function that simulates observations, given the animal's [`State`](@ref) and a `ModelObs` instance. This accepts the following arguments:
-   `state`: A [`State`](@ref) instance;
-   `model_obs`: A [`ModelObs`](@ref) instance;
-   `t`: An integer that defines the time step;

Methods are implemented for all built-in sub-types. Methods can be defined for new sub-types, such as `ModelObsDepthNormal`, as follows:
```
function Patter.simulate_obs(state::StateXYZD, model_obs::ModelObsDepthNormal, t::Int64)
    dbn   = truncated(Normal(state.z, model_obs.depth_sigma), 0, state.map_value)
    rand(dbn)
end
```

[`Patter.simulate_obs()`](@ref) is wrapped by [`simulate_yobs()`](@ref) for the simulation of observations.

# Log probabilities 

[`Patter.logpdf_obs()`](@ref) is a generic function that calculates the log probability (density) of an observation, given the animal's [`State`](@ref) and a `ModelObs` instance. This accepts the following arguments:
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

[`Patter.logpdf_obs()`](@ref) is used in [`particle_filter()`](@ref) to evaluate the log-probability of the data given particle samples.
        
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
@doc (@doc ModelObs) ModelObsAcousticLogisTrunc

function simulate_obs(state::State, model_obs::ModelObsAcousticLogisTrunc, t::Int64)
    # Evaluate the distance between the particle and receiver
    dist = distance(state.x, state.y, model_obs.receiver_x, model_obs.receiver_y)
    # Define probability of detection
    prob = ifelse(dist > model_obs.receiver_gamma, 0.0, logistic(model_obs.receiver_alpha + model_obs.receiver_beta * dist))
    # Define distribution
    rand(Bernoulli(prob)) + 0
end
@doc (@doc ModelObs) simulate_obs

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
@doc (@doc ModelObs) logpdf_obs


#########################
#########################
#### Uniform depth model 

struct ModelObsDepthUniform <: ModelObs
    sensor_id::Int64
    depth_shallow_eps::Float64
    depth_deep_eps::Float64
end
@doc (@doc ModelObs) ModelObsDepthUniform

function simulate_obs(state::State, model_obs::ModelObsDepthUniform, t::Int64)
    if model_obs.depth_shallow_eps == model_obs.depth_deep_eps == 0.0
        return state.map_value
    end 
    a     = max(0, state.map_value - model_obs.depth_shallow_eps)
    b     = state.map_value + model_obs.depth_deep_eps
    dbn   = Uniform(a, b)
    rand(dbn)
end 

function logpdf_obs(state::State, model_obs::ModelObsDepthUniform, t::Int64, obs::Float64)
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

struct ModelObsDepthNormalTrunc <: ModelObs 
    sensor_id::Int64
    # Normal distribution variance
    depth_sigma::Float64
    # The individual may be `depth_deep_eps` deeper than the bathymetric depth
    depth_deep_eps::Float64
end 
@doc (@doc ModelObs) ModelObsDepthNormalTrunc

# Truncated normal depth model log probability 
# * The probability is highest if the individual is on the seabed
# * The individual can be up to model_obs.depth_deep_eps deeper than the seabed
# * Probability decays away from the seabed toward the surface
function simulate_obs(state::State, model_obs::ModelObsDepthNormalTrunc, t::Int64)
    dbn   = truncated(Normal(state.map_value, model_obs.depth_sigma), 
                      0.0, state.map_value + model_obs.depth_deep_eps)
    rand(dbn)
end 

function logpdf_obs(state::State, model_obs::ModelObsDepthNormalTrunc, t::Int64, obs::Float64)
    dbn   = truncated(Normal(state.map_value, model_obs.depth_sigma), 
                      0.0, state.map_value + model_obs.depth_deep_eps)
    logpdf(dbn, obs)
end 