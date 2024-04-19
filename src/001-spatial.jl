import GeoArrays

export extract

# Extract values from an environmental layer
# * The default method assumes env is a GeoArray
# * A custom method is required to handle other types (e.g., shapefiles)

function extract(env::GeoArrays.GeoArray, x::Real, y::Real)
    cell = GeoArrays.indices(env, (x, y))
    if checkbounds(Bool, env, cell)
        env[cell]
    else
        eltype(env)(NaN)
    end
end

# Determine if a location is valid (i.e., in water)

# 2D case
function is_valid(env, x::Real, y::Real)
    z_env = extract(env, x, y)
    if isnan(z_env) 
        false
    else
        true
    end 
end

# 3D case 
function is_valid(env, x::Real, y::Real, z::Real)
    z_env = extract(env, x, y)
    if isnan(z_env) 
        false
    elseif 0 < z <= z_env
        true 
    else 
        false
    end 
end

# Distances

# 2D case
function distance(x0::Real, y0::Real, x1::Real, y1::Real)
    # sqrt((x0 - x1)^2 + (y0 - y1)^2)
    hypot(x0 - x1, y0 - y1)
end