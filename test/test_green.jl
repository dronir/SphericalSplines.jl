
G = SphericalSplines.Green

@testset "Green function" begin
    v = [1.0, 0.0, 0.0]
    @test G(v, v) ≈ 1.0 / 4π
    @test G(v, -v) ≈ 1.0 / 4π - π/24
    @test G(-v, v) ≈ 1.0 / 4π - π/24
    @test G(-v, -v) ≈ 1.0 / 4π 
end
