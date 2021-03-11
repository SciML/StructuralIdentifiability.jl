@testset "Check identifiability from given IO-equations" begin
    R, (x, y, z, t) = PolynomialRing(Nemo.QQ, ["x", "y", "z", "t"])

    @test check_identifiability(
        x * (y + z) * t^4 - x^3 * y * z * t^5 + 7 + t * x,
        [x, y, z],
        0.99
    ) == [true, false, false]

    @test check_identifiability(
        [z + t * (x + y) + z^3, z - t^2 * x * y *5 + t^2],
        [x, y],
        [x, y, x + y, x - y, x^2 + y^2],
        0.99
    ) == [false, false, true, false, true]

end
