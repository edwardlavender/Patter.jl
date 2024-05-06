using Patter
using Test

@testset "Patter tests" begin

    @testset "test-particle-filter.jl" begin
        include("test-particle-filter.jl")
    end

end