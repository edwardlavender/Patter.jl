using Distributions
using LogExpFunctions: logistic, log1pexp

export ModelObs
export ModelObsAcousticLogisTrunc, ModelObsAcousticLogisStd
export ModelObsDepthUniform, ModelObsDepthNormalTrunc, ModelObsDepthNormal


"""
    Observation models

# Structures

`ModelObs` is an abstract type used to hold parameters for observation models. 

`ModelObsAcousticLogisTrunc` is a `ModelObs` structure for an acoustic observation and a truncated logistic detection probability model:
        
- `x`, `y`: The coordinates of the receiver;
- `alpha`, `beta`, `gamma`: The parameters of a logistic detection probability model;

`ModelObsDepthUniform` is `ModelObs` structure for a depth observation and a uniform depth model, in which we assume the individual must be located in an envelope around the bathymetry depth, defined by two error terms:

- `depth_shallow_eps`
- `depth_deep_eps`

`ModelObsDepthNormalTrunc` is a `ModelObs` structure for a depth observation and a truncated normal model, in which we assume the likelihood of a depth observation given by a normal distribution centred at the bathymetric depth in a location and defined by the parameters:

- `sigma`: The standard deviation of the normal distribution;
- `deep_depth_eps`: The deep truncation parameter;

# Simulation

`simulate_obs()` is a generic function that simulates observations from a `ModelObs` instance. 

# Density 

`logpdf_obs()` is a generic function that calculates the log probability (density) of an observation, given the animal's `state` and a `ModelObs` instance.
        
"""
abstract type ModelObs end


#########################
#########################
#### Truncated logistic acoustic observation model

struct ModelObsAcousticLogisTrunc <: ModelObs
    # Receiver coordinates
    x::Float64
    y::Float64
    # Detection probability parameters
    alpha::Float64
    beta::Float64    
    gamma::Float64
end
@doc (@doc ModelObs) ModelObsAcousticLogisTrunc

function simulate_obs(state::State, model::ModelObsAcousticLogisTrunc, t::Int64)
    # Evaluate the distance between the particle and receiver
    dist = distance(state.x, state.y, model.x, model.y)
    # Define probability of detection
    prob = ifelse(dist > model.gamma, 0.0, logistic(model.alpha + model.beta * dist))
    # Define distribution
    rand(Bernoulli(prob))
end
@doc (@doc ModelObs) simulate_obs

function logpdf_obs(state::State, model::ModelObsAcousticLogisTrunc, t::Int64, obs::Int64)

    # Evaluate the distance between the particle and receiver
    dist = distance(state.x, state.y, model.x, model.y)

    # Calculate log probability given detection (1)
    if obs == 1
        if dist > model.gamma 
            return -Inf
        else 
            return -log1pexp(-(model.alpha + model.beta * dist))
        end 
    
    # Calculate log probability given non detection (0)
    elseif obs == 0
        if (dist > model.gamma)
            return 0.0
        else 
            return -log1pexp(model.alpha + model.beta * dist)
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
    depth_shallow_eps::Float64
    depth_deep_eps::Float64
end
@doc (@doc ModelObs) ModelObsDepthUniform

function simulate_obs(state::State, model::ModelObsDepthUniform, t::Int64)
    a     = max(0, state.map_value - model.depth_shallow_eps)
    b     = state.map_value + model.depth_deep_eps
    dbn   = Uniform(a, b)
    rand(dbn)
end 

function logpdf_obs(state::State, model::ModelObsDepthUniform, t::Int64, obs::Float64)
    a     = max(0, state.map_value - model.depth_shallow_eps)
    b     = state.map_value + model.depth_deep_eps
    dbn   = Uniform(a, b)
    logpdf(dbn, obs)
end 


#########################
#########################
#### Truncated normal depth model

struct ModelObsDepthNormalTrunc <: ModelObs 
    # Normal distribution variance
    sigma::Float64
    # The individual may be `deep_depth_eps` deeper than the bathymetric depth
    deep_depth_eps::Float64
end 
@doc (@doc ModelObs) ModelObsDepthNormalTrunc

# Truncated normal depth model log probability 
# * The probability is highest if the individual is on the seabed
# * The individual can be up to model.deep_depth_eps deeper than the seabed
# * Probability decays away from the seabed toward the surface
function simulate_obs(state::State, model::ModelObsDepthNormalTrunc, t::Int64)
    dbn   = truncated(Normal(state.map_value, model.sigma), 
                      0.0, state.map_value + model.deep_depth_eps)
    rand(dbn)
end 

function logpdf_obs(state::State, model::ModelObsDepthNormalTrunc, t::Int64, obs::Float64)
    dbn   = truncated(Normal(state.map_value, model.sigma), 
                      0.0, state.map_value + model.deep_depth_eps)
    logpdf(dbn, obs)
end 