@testset "Linear relations over the rationals" begin
    R, (a, b, c) = QQ["a", "b", "c"]

    cases = []

    push!(
        cases,
        Dict(
            :rff => StructuralIdentifiability.RationalFunctionField([a + 9]),
            :correct => Set([a // one(R)]),
            :degree => 2,
        ),
    )

    push!(
        cases,
        Dict(
            :rff => StructuralIdentifiability.RationalFunctionField([
                a * b // R(1),
                (b * c + a * b) // (a * b),
            ]),
            :correct => Set([a * b // one(R), b * c // one(R)]),
            :degree => 2,
        ),
    )

    f = [a^2 + b^2, a^3 + b^3, a^4 + b^4]
    rff = StructuralIdentifiability.RationalFunctionField(f)
    push!(cases, Dict(:rff => rff, :correct => Set([a + b // one(R)]), :degree => 1))
    push!(
        cases,
        Dict(
            :rff => rff,
            :correct => Set([a + b // one(R), a * b // one(R)]),
            :degree => 2,
        ),
    )

    f = [9a^7 + 10b^6, b^10 - 5b^2]
    rff = StructuralIdentifiability.RationalFunctionField(f)
    push!(cases, Dict(:rff => rff, :correct => Set(empty([a // one(R)])), :degree => 1))
    push!(
        cases,
        Dict(:rff => rff, :correct => Set([a^7 + (10 // 9) * b^6 // one(R)]), :degree => 7),
    )
    push!(
        cases,
        Dict(
            :rff => rff,
            :correct => Set([a^7 + (10 // 9) * b^6 // one(R), b^10 - 5b^2 // one(R)]),
            :degree => 12,
        ),
    )

    push!(
        cases,
        Dict(
            :rff => StructuralIdentifiability.RationalFunctionField([a, a * b + b * c]),
            :correct => Set([a // one(R), a * b + b * c // one(R)]),
            :degree => 2,
        ),
    )

    for c in cases
        @test Set(StructuralIdentifiability.polynomial_generators(c[:rff], c[:degree])) ==
              c[:correct]
    end

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
    relations = StructuralIdentifiability.polynomial_generators(rff, 2)
    @test (x1 * p4 + p2 * x2) // one(R) in relations

    ###
    # Some arbitrary generators for the SLIQR model
    R, (b, e, In, S, Ninv, s, Q, g, u, a, y, L) =
        polynomial_ring(QQ, [:b, :e, :In, :S, :Ninv, :s, :Q, :g, :u, :a, :y, :L])
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
    relations = StructuralIdentifiability.polynomial_generators(rff, 2)
    @test s * Q - Q * a in relations
end
