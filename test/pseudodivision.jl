@testset "Pseudodivision" begin
    R, (x, y, z) = Nemo.polynomial_ring(Nemo.QQ, ["x", "y", "z"])

    @test iszero(pseudodivision(x^2 - y^2, x - y, x))
    @test pseudodivision(y * x^4 + x^3 + x, z * x^2, x) == x
end
