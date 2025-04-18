# Run tests for selected Patter.jl functions
# (Many functions are also tested via the patter R package)

using Patter
using Test

@testset "Patter.jl" begin

    include("test-001-utilities.jl")
    include("test-002-spatial.jl")
    include("test-009-particle-filter.jl")
    include("test-012-Julia-from-R.jl")

end