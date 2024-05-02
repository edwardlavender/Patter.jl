using Random

"""
    # Systematic resampling algorithm
    
    Given the weight vector w, resample a set of *indices* based on low-variance resampling algorithm from Thrun, Burgard, and Fox's "Probabilistic Robotics".
    
    # Source 

    Code adapted from https://github.com/JuliaStats/StatsBase.jl/issues/124.
    
    # Example
    
    ```Julia
    X = ["A", "B", "C", "D"]
    w = [0, 0, 0.75, 0.25]

    idx = resample(w, 12)
    X[idx]
    ```
    """
function resample(w::Vector{Float64}, n::Int = length(w))
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