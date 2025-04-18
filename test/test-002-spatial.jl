@testset "Patter.is_valid()" begin

    @test Patter.is_valid(100)
    @test !Patter.is_valid(NaN)

    @test Patter.is_valid(101, 100)
    @test !Patter.is_valid(NaN, 100)
    @test !Patter.is_valid(10, 100)
    @test !Patter.is_valid(100, 101)
    
end 