@testset "Patter.resample()" begin

    w = [20.0, 10 ,80]
    ii = Patter.resample(w, Int(sum(w)))

    @test sum(ii .== 1) == w[1]
    @test sum(ii .== 2) == w[2]
    @test sum(ii .== 3) == w[3]

    for n in 0:10
        @test length(Patter.resample(w, n)) == n
    end

end