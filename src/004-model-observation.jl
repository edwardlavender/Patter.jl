using Distributions
using LogExpFunctions: logistic, log1pexp

export ModelObs
export ModelObsAcousticLogisTrunc, ModelObsAcousticLogisStd
export ModelObsDepthUniform, ModelObsDepthNormalTrunc, ModelObsDepthNormal


"""
    Observation models

# Structures

`ModelObs` is an abstract type used to hold parameters for observation models. For communication with `R`, all subtypes should include a `sensor_id` field. 

`ModelObsAcousticLogisTrunc` is a `ModelObs` structure for an acoustic observation and a truncated logistic detection probability model:
 
- `sensor_id`: The receiver ID;
- `receiver_x`, `receiver_y`: The coordinates of the receiver;
- `receiver_alpha`, `receiver_beta`, `receiver_gamma`: The parameters of a logistic detection probability model;

`ModelObsDepthUniform` is `ModelObs` structure for a depth observation and a uniform depth model, in which we assume the individual must be located in an envelope around the bathymetry depth, defined by two error terms:

- `sensor_id`: The sensor ID;
- `depth_shallow_eps`
- `depth_deep_eps`

`ModelObsDepthNormalTrunc` is a `ModelObs` structure for a depth observation and a truncated normal model, in which we assume the likelihood of a depth observation given by a normal distribution centred at the bathymetric depth in a location and defined by the parameters:

- `sensor_id`: The sensor ID;
- `depth_sigma`: The standard deviation of the normal distribution;
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

function simulate_obs(state::State, model::ModelObsAcousticLogisTrunc, t::Int64)
    # Evaluate the distance between the particle and receiver
    dist = distance(state.x, state.y, model.receiver_x, model.receiver_y)
    # Define probability of detection
    prob = ifelse(dist > model.receiver_gamma, 0.0, logistic(model.receiver_alpha + model.receiver_beta * dist))
    # Define distribution
    rand(Bernoulli(prob)) + 0
end
@doc (@doc ModelObs) simulate_obs

function logpdf_obs(state::State, model::ModelObsAcousticLogisTrunc, t::Int64, obs::Int64)

    # Evaluate the distance between the particle and receiver
    dist = distance(state.x, state.y, model.receiver_x, model.receiver_y)

    # Calculate log probability given detection (1)
    if obs == 1
        if dist > model.receiver_gamma 
            return -Inf
        else 
            return -log1pexp(-(model.receiver_alpha + model.receiver_beta * dist))
        end 
    
    # Calculate log probability given non detection (0)
    elseif obs == 0
        if dist > model.receiver_gamma
            return 0.0
        else 
            return -log1pexp(model.receiver_alpha + model.receiver_beta * dist)
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
    sensor_id::Int64
    # Normal distribution variance
    depth_sigma::Float64
    # The individual may be `depth_deep_eps` deeper than the bathymetric depth
    depth_deep_eps::Float64
end 
@doc (@doc ModelObs) ModelObsDepthNormalTrunc

# Truncated normal depth model log probability 
# * The probability is highest if the individual is on the seabed
# * The individual can be up to model.depth_deep_eps deeper than the seabed
# * Probability decays away from the seabed toward the surface
function simulate_obs(state::State, model::ModelObsDepthNormalTrunc, t::Int64)
    dbn   = truncated(Normal(state.map_value, model.depth_sigma), 
                      0.0, state.map_value + model.depth_deep_eps)
    rand(dbn)
end 

function logpdf_obs(state::State, model::ModelObsDepthNormalTrunc, t::Int64, obs::Float64)
    dbn   = truncated(Normal(state.map_value, model.depth_sigma), 
                      0.0, state.map_value + model.depth_deep_eps)
    logpdf(dbn, obs)
end 