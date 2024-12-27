using DataFrames

@testset "Patter.julia_get_*()" begin
    
    # Patter.julia_get_xinit()
    d = DataFrame(map_value = [0, 1], x = [1, 2], y = [3, 4])
    @test Patter.julia_get_xinit(StateXY, d) == [StateXY(0, 1, 3), StateXY(1, 2, 4)]

    # Patter.julia_get_model_obs_types()
    @test Patter.julia_get_model_obs_types(["ModelObsAcousticLogisTrunc"]) == [ModelObsAcousticLogisTrunc]
    @test Patter.julia_get_model_obs_types(["ModelObsAcousticLogisTrunc", "ModelObsDepthUniform"]) == [ModelObsAcousticLogisTrunc, ModelObsDepthUniform]

    # Patter.r_get_states()
    # * Define example State matrix ()
    # - Two rows: two particles
    # - Three columns: three time steps 
    timesteps  = [1, 2, 3]
    timestamps = [DateTime("2016-01-01T00:00:00") + Minute(2) * i for i in 0:2]
    state      = [StateXY(0.0, 1.0, 2.0)  StateXY(0.0, 3.0, 4.0)  StateXY(0.0, 5.0, 6.0);
                  StateXY(0.0, 7.0, 8.0) StateXY(0.0, 9.0, 10.0)    StateXY(0.0, 11.0, 12.0)]
    # * Convert to DataFrame
    output    = Patter.r_get_states(state, timesteps, timestamps)
    # * Define epected output
    expected = DataFrame(
        path_id = [1, 1, 1, 2, 2, 2],
        timestep = [1, 2, 3, 1, 2, 3],
        timestamp = [
            DateTime("2016-01-01T00:00:00"),
            DateTime("2016-01-01T00:02:00"),
            DateTime("2016-01-01T00:04:00"),
            DateTime("2016-01-01T00:00:00"),
            DateTime("2016-01-01T00:02:00"),
            DateTime("2016-01-01T00:04:00")
        ],
        map_value = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        x = [1.0, 3.0, 5.0, 7.0, 9.0, 11.0],
        y = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
    )
    # * Test
    @test output == expected

    # Patter.r_get_states(): Abig matrix of StateXY objects
    # np = 1000
    # nt = 20000
    # state = [StateXY(rand(), rand(), rand()) for _ in 1:np, _ in 1:nt]
    # r_get_states(state)

end

