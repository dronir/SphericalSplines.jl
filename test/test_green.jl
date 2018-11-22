
G = SphericalSplines.Green

@testset "Green function" begin
    v = [1.0, 0.0, 0.0]
    @test G(v, v) ≈ 1.0 / 4π
    @test G(v, -v) ≈ 1.0 / 4π - π/24
    @test G(-v, v) ≈ 1.0 / 4π - π/24
    @test G(-v, -v) ≈ 1.0 / 4π 
    @test isapprox(G(v, [0.0, 1.0, 0.0]), (1 - π^2/12 + log(2)^2/2 - log(2)^2) / 4π; atol=TOL_DIGITS)
end
