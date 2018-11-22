
Li2 = SphericalSplines.Li2

@testset "Dilogarithm approximation to $TOL_DIGITS digits" begin
    @test isapprox(Li2(0.0), 0.0 ; atol=TOL_DIGITS)
    @test isapprox(Li2(-1.0), -π^2 / 12 ; atol=TOL_DIGITS)
    @test isapprox(Li2(0.5), π^2/12 - log(2)^2/2; atol=TOL_DIGITS)
    @test isapprox(Li2(1.0), π^2 / 6; atol=TOL_DIGITS)
end
