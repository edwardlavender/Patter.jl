using Distributions
using LogExpFunctions: log1pexp

export ModelObs
export ModelObsAcousticLogisTrunc, ModelObsAcousticLogisStd
export ModelObsDepthUniform, ModelObsDepthNormalTrunc, ModelObsDepthNormal

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

function log_prob_obs(state::State, model::ModelObsAcousticLogisTrunc, t::Int64, obs::Int64)

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


#########################
#########################
#### Standard logistic acoustic observation model parameters

struct ModelObsAcousticLogisStd
       # Receiver coordinates
       x::Float64
       y::Float64
       # Detection probability parameters
       d50::Float64
       k::Float64    
end 

function log_prob_obs(state::State, model::ModelObsAcousticLogisStd, t::Int64, obs::Int64)

    # Evaluate the distance between the particle and receiver
    dist = distance(state.x, state.y, model.x, model.y)

    # Calculate log probability given detection (1)
    if obs == 1
        return -log1pexp(-(model.d50 - dist) / model.k)
    
    # Calculate log probability given non detection (0)
    elseif obs == 0
        if dist > 15 * model.k + model.d50
            return zero(Float64)
        else 
            return -log1pexp((model.d50 - dist) / model.k)
        end

    else 
      error("Acoustic observations should be coded as 0 or 1.")  
    end

end


#########################
#########################
#### Uniform depth model 

# TO DO
# * Review the inclusion of the env dataset here as well

struct ModelObsDepthUniform{T} <: ModelObs
    env::T
    depth_shallow_eps::Float64
    depth_deep_eps::Float64
end

function log_prob_obs(state::State, model::ModelObsDepthUniform, t::Int64, obs::Float64)
    z_env = extract(model.env, state.x, state.y)
    if isnan(z_env) 
        return -Inf 
    end 
    a = max(0, z_env - model.depth_shallow_eps)
    b = z_env + model.depth_deep_eps
    logpdf(Uniform(a, b), obs)
end 


#########################
#########################
#### Truncated normal depth model

struct ModelObsDepthNormalTrunc{T} <: ModelObs 
    env::T
    # Normal distribution variance
    sigma::Float64
    # The individual may be `deep_depth_eps` deeper than the bathymetric depth
    deep_depth_eps::Float64
end 

# Truncated normal depth model log probability 
# * The probability is highest if the individual is on the seabed
# * The individual can be up to model.deep_depth_eps deeper than the seabed
# * Probability decays away from the seabed toward the surface
function log_prob_obs(state::State, model::ModelObsDepthNormalTrunc, t::Int64,  obs::Float64)
    z_env = extract(model.env, state.x, state.y)
    if isnan(z_env) 
        return -Inf 
    end 
    logpdf(truncated(Normal(z_env, model.sigma), 
                     0.0, z_env + model.deep_depth_eps), obs)
end 


#########################
#########################
#### Normal depth model

# Normal depth model 
struct ModelObsDepthNormal <: ModelObs
    # Normal distribution variance
    sigma::Float64
end

# Normal depth model log probability 
function log_prob_obs(state::StateXYZD, model::ModelObsDepthNormal, t::Int64, obs::Float64)
    logpdf(Normal(state.z, model.sigma), obs)
end 