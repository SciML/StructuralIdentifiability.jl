if GROUP == "All" || GROUP == "Core"
    @testset "RationalFunctionField" begin
        p = 0.99
        R, (a, b, c) = QQ["a", "b", "c"]

        f1 = [a, b, a + b + c]
        f2 = [2c, 3b, 5a]
        rff1 = StructuralIdentifiability.RationalFunctionField(f1)
        rff2 = StructuralIdentifiability.RationalFunctionField(f2)
        @test StructuralIdentifiability.fields_equal(rff1, rff2, p)
        @test all(StructuralIdentifiability.field_contains(rff1, [zero(R), one(R)], p))
        @test all(StructuralIdentifiability.field_contains(rff1, [7a, 9b, 11c], p))
        @test all(StructuralIdentifiability.field_contains(rff2, [7a, 9b, 11c], p))

        s1, s2 = a + b + c, a^2 + b^2 + c^2
        f1 = [s1, s2]
        rff1 = StructuralIdentifiability.RationalFunctionField(f1)
        @test !any(
            StructuralIdentifiability.field_contains(rff1, [a, b + c, a * b + b * c], p),
        )
        @test all(
            StructuralIdentifiability.field_contains(rff1, [a * b + b * c + a * c], p),
        )
        @test all(StructuralIdentifiability.field_contains(rff1, [(s1)^8 - (s2)^9 + 89], p))

        # Example in Section 5 from
        # https://mediatum.ub.tum.de/doc/685465/685465.pdf
        R, (x1, x2) = QQ["x1", "x2"]
        g1 = (x1^3 + x1 * x2 - 2) // (x1^2 - x2 - 1)
        g2 = (x1^2 + x1^2 * x2 + 7) // (x1 - x1^2 * x2^2)
        g3 = x1^2 + 3x1 * x2
        g4 = x1 * x2^2 + 5x1 * x2
        g5 = x1^3 * x2 - x2
        rff1 = StructuralIdentifiability.RationalFunctionField([g1, g2, g3, g4, g5])
        rff2 = StructuralIdentifiability.RationalFunctionField([x1, x2])
        @test StructuralIdentifiability.fields_equal(rff1, rff2, p)
    end
end
