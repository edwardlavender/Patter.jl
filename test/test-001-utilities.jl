using Dates
using DataFrames

@testset "Patter.check_names()" begin
    @test isnothing(Patter.check_names(DataFrame(x = 1, y = 2), "x"))
    @test_throws ErrorException Patter.check_names(DataFrame(x = 1, y = 2), "z")
    @test_throws ErrorException Patter.check_names(DataFrame(x = 1, y = 2), ["a", "b"])
end

@testset "Patter.diffsecs()" begin
    t1 = DateTime("2024-09-01T12:00:00")
    t2 = DateTime("2024-09-01T12:05:00")
    @test Patter.diffsecs(t2, t1) == 5 * 60.0

    t1 = DateTime("2024-09-01T12:00:00")
    t2 = DateTime("2024-09-01T11:55:00")
    @test Patter.diffsecs(t2, t1) == -5 * 60.0

    t1 = DateTime("2024-09-01T12:00:00")
    t2 = DataFrame(timeline = [DateTime("2024-09-01T12:10:00"), DateTime("2024-09-01T12:15:00")])
    @test Patter.diffsecs(t2.timeline, t1) == [10 * 60.0, 15 * 60.0]
end