using DataFrames
using Tables

"""
# Julia from R

A collection of functions that facilitate the translation of inputs from R into Julia.

# Details

* [`julia_get_xinit()`](@ref) gets a Vector of initial States from a DataFrame;
* [`julia_get_model_types()`](@ref) gets a Vector of `ModelObs` subtypes from a Vector of Strings;
* [`julia_get_models()`](@ref) gets a Vector of model instances from a Vector of DataFrames that contain parameters and a corresponding vector of model types;

"""
function julia_get end 

#### Get the Vector of initial States 
function julia_get_xinit(state::Type{<:State}, d::DataFrame)
    fields = fieldnames(state)
    [state((d[i, Symbol(field)] for field in fields)...) for i in 1:nrow(d)]
end 

# d = DataFrame(map_value = [0, 1], x = [1, 2], y = [3, 4])
# julia_get_xinit(StateXY, d)

#### Get a vector of ModelObs types 

function julia_get_model_types(models::Union{String, Vector{String}})
    models = isa(models, String) ? [models] : models
    return [getfield(Main, Symbol(model)) for model in models]
end

# julia_get_models(["ModelObsAcousticLogisTrunc"])
# julia_get_models(["ModelObsAcousticLogisTrunc", "ModelObsDepthUniform"])

#### Get a vector of ModelObs instances
function julia_get_models(parameters::Vector, model_types::Vector{DataType})
    sensors = Vector{ModelObs}()
    for (dataset, model_type) in zip(parameters, model_types)
        # TO DO
        # Force correct column order with a warning
        for row in Tables.namedtupleiterator(dataset)
            push!(sensors, model_type(row...))
        end 
    end 
    sensors
end 