using DataFrames
using Tables

# """
# # `Julia` from `R`

# A collection of internal functions that facilitate the translation of objects from `R` into `Julia` objects.

# # Details

# * [`Patter.julia_get_xinit()`](@ref) gets a `Vector` of initial [`State`](@ref)s from a `DataFrame`;
# * [`Patter.julia_get_model_obs_types()`](@ref) gets a `Vector` of [`ModelObs`](@ref) sub-types from a Vector of `String`s;
# * [`Patter.julia_get_model_obs()`](@ref) gets a `Vector` of [`ModelObs`](@ref) instances from a `Vector` of `DataFrame`s that contain parameters and a corresponding vector of [`ModelObs`](@ref) sub-types;

# """
# function julia_get end 

# Get the Vector of initial States 
function julia_get_xinit(state::Type{<:State}, d::DataFrame)
    fields = fieldnames(state)
    [state((d[i, Symbol(field)] for field in fields)...) for i in 1:nrow(d)]
end 

# Get a vector of ModelObs types 
function julia_get_model_obs_types(models::Union{String, Vector{String}})
    models = isa(models, String) ? [models] : models
    return [getfield(Main, Symbol(model)) for model in models]
end

# Get a vector of ModelObs instances
function julia_get_model_obs(parameters::Vector, model_types::Vector{DataType})
    sensors = Vector{ModelObs}()
    for (dataset, model_type) in zip(parameters, model_types)
        # TO DO
        # Force correct column order with a warning & check names
        for row in Tables.namedtupleiterator(dataset)
            push!(sensors, model_type(row...))
        end 
    end 
    sensors
end 