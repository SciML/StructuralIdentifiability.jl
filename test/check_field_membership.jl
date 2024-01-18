if GROUP == "All" || GROUP == "Core"
    @testset "Check field membership" begin
        R, (x, y, z) = Nemo.PolynomialRing(Nemo.QQ, ["x", "y", "z"])

        @test field_contains(
            RationalFunctionField([[R(1), x + y], [R(1), x * y], [z, (x + y)^2]]),
            [(x^2 + y^2) // R(1), (x^3 + y^3) // (z - x * y), R(1) // (z + x + y), z // x],
            0.99,
        ) == [true, true, true, false]

        @test field_contains(
            RationalFunctionField([[
                x + y + z,
                x^2 + y^2 + z^2,
                (x + y + z)^2,
                x^3 + y^3 + z^3,
            ]]),
            [x + y + z // 1, x * y * z // 1, x + y + 2 * z // 1, x // (y + z)],
            0.99,
        ) == [true, true, false, false]

        @test field_contains(
            RationalFunctionField([
                x + y + z // 1,
                x * y + y * z + z * x // 1,
                x * y * z // 1,
            ]),
            [x^2 + y^2 + z^2, x^6 + y^6 + z^6, x - y + z, x^2 - y^2 + z^2],
            0.99,
        ) == [true, true, false, false]
    end
end
