using DataFrames
using Tables

"""
# Julia from R

A collection of functions that facilitate the translation of inputs from R into Julia.

# Details

* [`julia_get_xinit()`] gets a Vector of initial States from a DataFrame;
* [`julia_get_model_types()`] gets a Vector of `ModelObs` subtypes from a Vector of Strings;
* [`julia_get_models()`] gets a Vector of model instances from a Vector of DataFrames that contain parameters and a corresponding vector of model types;

"""
function julia_get end 

#### Get the Vector of initial States 

function julia_get_xinit(state::Type{StateXY}, d::DataFrame)
    [StateXY(d.map_value[i], d.x[i], d.y[i]) for i in 1:nrow(d)]
end
@doc (@doc julia_get) julia_get_xinit

function julia_get_xinit(state::Type{StateXYZ}, d::DataFrame)
    [StateXYZ(d.map_value[i], d.x[i], d.y[i], d.z[i]) for i in 1:nrow(d)]
end

function julia_get_xinit(state::Type{StateXYZD}, d::DataFrame)
    [StateXYZD(d.map_value[i], d.x[i], d.y[i], d.z[i], d.angle[i]) for i in 1:nrow(d)]
end

# d = DataFrame(map_value = [0, 1], x = [1, 2], y = [3, 4])
# julia_get_xinit(StateXY, d)

#### Get a vector of ModelObs types 

function julia_get_model_types(models::Vector{String})
    [eval(Symbol(model)) for model in models]
end

# julia_get_models(["ModelObsAcousticLogisTrunc"])
# julia_get_models(["ModelObsAcousticLogisTrunc", "ModelObsDepthUniform"])

#### Get a vector of ModelObs instances
function julia_get_models(parameters::Vector, model_types::Vector{DataType})
    sensors = Vector{ModelObs}()
    for (dataset, model_type) in zip(parameters, model_types)
        for row in Tables.namedtupleiterator(dataset)
            push!(sensors, model_type(row...))
        end 
    end 
    sensors
end 