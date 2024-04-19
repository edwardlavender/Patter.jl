import Random

export pf_init

# Low variance resampling algorithm 
function pf_resample_lv(w::Vector{Float64}, n::Int = length(w))
    w = w ./ sum(w)
    idx = zeros(Int, n)
    r = rand() * 1/n
    c = w[1]
    i = 1
    for m in 1:n
        # Calculate the next sample point
        U = r + (m - 1) * (1/n)
        # Find the first weight that puts us past the sample point
        while c < U
            i += 1
            c = c + w[i]
        end
        idx[m] = i
    end
    Random.shuffle!(idx)
    idx
end

# Initial particles 
function pf_init(state::Vector, model::ModelMove, n::Int, xlim = nothing, ylim = nothing)

    # (optional) Define xlim and ylim
    env = model.env
    bb = GeoArrays.bbox(env)
    if isnothing(xlim)
        xlim = (bb.min_x, bb.max_x)
    end
    if isnothing(ylim)
        ylim = (bb.min_y, bb.max_y)
    end
    
    # Sample within limits
    # * We assume that there are at least some valid locations within xlim & ylim
    xinit = state
    while length(xinit) < n
        pinit = rinit(state, model, xlim, ylim)
        valid = is_valid(model.env, pinit.x, pinit.y)
        if valid
            push!(xinit, pinit)
        end 
    end
    
    xinit

    end