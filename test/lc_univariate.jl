@testset "Univariate leading coefficient" begin
    R, (x, y, z) = Nemo.polynomial_ring(Nemo.QQ, ["x", "y", "z"])
    p = x^2 * y + x^2 * (z + z^3) + y - 5
    @test lc_univariate(p, x) == y + z + z^3
    @test lc_univariate(p, z) == x^2
    @test lc_univariate(p, y) == 1 + x^2
    @test lc_univariate(x + y, z) == x + y
end
