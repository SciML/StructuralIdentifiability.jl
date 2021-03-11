@testset "Check field membership" begin
    R, (x, y, z) = PolynomialRing(Nemo.QQ, ["x", "y", "z"])

    @test check_field_membership(
        [[R(1), x + y], [R(1), x * y], [z, (x + y)^2]],
        [(x^2 + y^2) // R(1), (x^3 + y^3) // (z - x * y), R(1) // (z + x + y), z // x],
        0.99
    ) == [true, true, true, false]

    @test check_field_membership(
        [[x + y + z, x^2 + y^2 + z^2, (x + y + z)^2, x^3 + y^3 + z^3]],
        [x + y + z, x * y * z, x + y + 2 * z, x // (y + z)],
        0.99
    ) == [true, true, false, false]
end
