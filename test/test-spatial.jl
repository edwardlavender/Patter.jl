@testset "is_valid" begin

    @test !is_valid(NaN)
    @test is_valid(100)
    
    @test !is_valid(NaN, 100)
    @test is_valid(10, 100)
    @test !is_valid(101, 100)
    
end 