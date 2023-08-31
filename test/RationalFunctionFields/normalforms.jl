
eq_up_to_the_order(a, b) = issubset(a, b) && issubset(b, a)

@testset "Linear relations over the rationals" begin
    for strategy in (:deterministic, :monte_carlo)
        R, (a, b, c) = QQ["a", "b", "c"]

        f = [a + 9]
        rff = StructuralIdentifiability.RationalFunctionField(f)
        relations = StructuralIdentifiability.monomial_generators_up_to_degree(
            rff,
            2,
            strategy = strategy,
        )
        @test eq_up_to_the_order(relations, [a, a^2])

        f = [a * b // R(1), (b * c + a * b) // (a * b)]
        rff = StructuralIdentifiability.RationalFunctionField(f)
        relations = StructuralIdentifiability.monomial_generators_up_to_degree(
            rff,
            2,
            strategy = strategy,
        )
        @test eq_up_to_the_order(relations, [a * b // R(1), b * c // R(1)])

        R, (a, b, c) = QQ["a", "b", "c"]
        f = [a^2 + b^2, a^3 + b^3, a^4 + b^4]
        rff = StructuralIdentifiability.RationalFunctionField(f)
        relations = StructuralIdentifiability.monomial_generators_up_to_degree(
            rff,
            1,
            strategy = strategy,
        )
        @test eq_up_to_the_order(relations, [a + b])
        relations = StructuralIdentifiability.monomial_generators_up_to_degree(
            rff,
            2,
            strategy = strategy,
        )
        @test eq_up_to_the_order(relations, [a + b, a * b, a^2 + b^2])

        strategy == :deterministic && continue
        f = [9a^7 + 10b^6, b^10 - 5b^2]
        rff = StructuralIdentifiability.RationalFunctionField(f)
        relations = StructuralIdentifiability.monomial_generators_up_to_degree(
            rff,
            1,
            strategy = strategy,
        )
        @test eq_up_to_the_order(relations, empty(f))
        relations = StructuralIdentifiability.monomial_generators_up_to_degree(
            rff,
            7,
            strategy = strategy,
        )
        @test eq_up_to_the_order(relations, [a^7 + (10 // 9) * b^6])
        relations = StructuralIdentifiability.monomial_generators_up_to_degree(
            rff,
            12,
            strategy = strategy,
        )
        @test eq_up_to_the_order(relations, [a^7 + (10 // 9) * b^6, b^10 - 5b^2])
    end

    # Regression test.
    # LV model.
    R, (x1, p2, p4, y1, x2, x3, u, p1, p3) =
        QQ["x1", "p2", "p4", "y1", "x2", "x3", "u", "p1", "p3"]
    f = [
        x3 // one(R),
        x2 * x1 // one(R),
        p1 * p3 // one(R),
        p2 * p4 // one(R),
        p1 + p3 // one(R),
        (p2 * x2 + p4 * x1) // (x2 * x1),
        (p2 * x2 - p4 * x1) // (p1 - p3),
    ]

    rff = StructuralIdentifiability.RationalFunctionField(f)
    relations = StructuralIdentifiability.monomial_generators_up_to_degree(
        rff,
        2,
        strategy = :monte_carlo,
    )
    @test (x1 * p4 + p2 * x2) // one(R) in relations
end
