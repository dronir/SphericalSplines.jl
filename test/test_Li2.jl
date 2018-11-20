
Li2 = SphericalSplines.Li2

const DIGITS=5

@testset "Dilogarithm approximation to $DIGITS digits" begin
    @test isapprox(Li2(0.0), 0.0 ; atol=DIGITS)
    @test isapprox(Li2(-1.0), -π^2 / 12 ; atol=DIGITS)
    @test isapprox(Li2(0.5), π^2/12 - log(2)^2/2; atol=DIGITS)
    @test isapprox(Li2(1.0), π^2 / 6; atol=DIGITS)
end
