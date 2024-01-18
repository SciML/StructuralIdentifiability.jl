if GROUP == "All" || GROUP == "Core"
    @testset "Coefficient extraction for rational fucntions" begin
        R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
        C = extract_coefficients_ratfunc(
            (x^2 + y * z - y^2 * z^3 + 3 * x * z^3) // (x + y + z + z^2 * (x^2 + 1)),
            [z],
        )

        @test Set(C) == Set([
            one(R) // 1,
            (3 * x - y^2) // 1,
            y // 1,
            x^2 // 1,
            (x + y) // 1,
            (x^2 + 1) // 1,
        ])

        R, (x, y) = PolynomialRing(QQ, ["x", "y"])
        f = (x^2 + y^2) // (1 - x - 3 * y)
        @test Set(extract_coefficients_ratfunc(f, Vector{Nemo.fmpq_mpoly}())) ==
              Set([f, one(R) // 1])

        R, (x, y, u, v) = PolynomialRing(QQ, ["x", "y", "u", "v"])
        C = extract_coefficients_ratfunc(
            (x + (y + 3) * u * v + y^2 * v^3) // (u + 3 * v - (x^2 + y^2) * u^2),
            [u, v],
        )
        @test Set(C) == Set([
            x // 1,
            (y + 3) // 1,
            y^2 // 1,
            one(R) // 1,
            3 * one(R) // 1,
            -(x^2 + y^2) // 1,
        ])
    end

    @testset "Coefficient extraction for polynomials" begin
        R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
        C = extract_coefficients((y + z + 8), [x])
        R_coef = parent(first(values(C)))
        y, z = gens(R_coef)
        @test symbols(R_coef) == [:y, :z]
        @test C == Dict([0] => y + z + 8)

        R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
        C = extract_coefficients((x^2 + y * z - y^2 * z^3 + 3 * x * z^3), [z])
        R_coef = parent(first(values(C)))
        x, y = gens(R_coef)
        @test symbols(R_coef) == [:x, :y]
        @test C == Dict([3] => 3x - y^2, [1] => y, [0] => x^2)

        R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
        C = extract_coefficients((x^2 + y * z - y^2 * z^3 + 3 * x * z^3), [x, z])
        R_coef = parent(first(values(C)))
        y = gens(R_coef)[1]
        @test symbols(R_coef) == [:y]
        @test C ==
              Dict([1, 3] => R_coef(3), [2, 0] => R_coef(1), [0, 3] => -y^2, [0, 1] => y)
    end
end
