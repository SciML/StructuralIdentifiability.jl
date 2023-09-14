eq_up_to_the_order(a, b) = issubset(a, b) && issubset(b, a)

@testset "Linear relations over the rationals" begin
    R, (a, b, c) = QQ["a", "b", "c"]

    f = [a + 9]
    rff = StructuralIdentifiability.RationalFunctionField(f)
    relations = StructuralIdentifiability.monomial_generators_up_to_degree(rff, 2)
    @test eq_up_to_the_order(relations, [a])

    f = [a * b // R(1), (b * c + a * b) // (a * b)]
    rff = StructuralIdentifiability.RationalFunctionField(f)
    relations = StructuralIdentifiability.monomial_generators_up_to_degree(rff, 2)
    @test eq_up_to_the_order(relations, [a * b // R(1), b * c // R(1)])

    R, (a, b, c) = QQ["a", "b", "c"]
    f = [a^2 + b^2, a^3 + b^3, a^4 + b^4]
    rff = StructuralIdentifiability.RationalFunctionField(f)
    relations = StructuralIdentifiability.monomial_generators_up_to_degree(rff, 1)
    @test eq_up_to_the_order(relations, [a + b])
    relations = StructuralIdentifiability.monomial_generators_up_to_degree(rff, 2)
    @test eq_up_to_the_order(relations, [a + b, a * b, a^2 + b^2])

    f = [9a^7 + 10b^6, b^10 - 5b^2]
    rff = StructuralIdentifiability.RationalFunctionField(f)
    relations = StructuralIdentifiability.monomial_generators_up_to_degree(rff, 1)
    @test eq_up_to_the_order(relations, empty(f))
    relations = StructuralIdentifiability.monomial_generators_up_to_degree(rff, 7)
    @test eq_up_to_the_order(relations, [a^7 + (10 // 9) * b^6])
    relations = StructuralIdentifiability.monomial_generators_up_to_degree(rff, 12)
    @test eq_up_to_the_order(relations, [a^7 + (10 // 9) * b^6, b^10 - 5b^2])

    # Regression tests
    ###
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

    ###
    R, (a, b, c) = QQ["a", "b", "c"]
    f = [a, a * b + b * c]
    rff = StructuralIdentifiability.RationalFunctionField(f)
    relations = StructuralIdentifiability.monomial_generators_up_to_degree(rff, 2)
    @test eq_up_to_the_order(relations, [a, a * b + b * c])

    ###
    # Some arbitrary generators for the SLIQR model
    R, (b, e, In, S, Ninv, s, Q, g, u, a, y, L) =
        PolynomialRing(QQ, [:b, :e, :In, :S, :Ninv, :s, :Q, :g, :u, :a, :y, :L])
    f = [
        In // one(R),
        s // one(R),
        Ninv // one(R),
        b // one(R),
        (g + a) // one(R),
        (e * s * g - s * g + g * a) // one(R),
        (e * S - S) // (e * Q),
        (e * S * s - S * s + S * a) // e,
        (s * Q^2 - Q^2 * a) // (e * g - g),
        (e * In + e * L - In - Q - L) // (e * Q),
    ]
    rff = StructuralIdentifiability.RationalFunctionField(f)
    relations = StructuralIdentifiability.monomial_generators_up_to_degree(rff, 2)
    @test s * Q - Q * a in relations
end

using Test, Logging, Nemo
Base.global_logger(ConsoleLogger(Logging.Info))

R, (a, b, c) = QQ["a", "b", "c"]

f = [a, a * b + b * c]
fracs = StructuralIdentifiability.beautifuly_generators(
    StructuralIdentifiability.RationalFunctionField(f),
)

relations = StructuralIdentifiability.linear_relations_between_normal_forms(fracs, 2)
