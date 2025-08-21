@testset "RationalFunctionField: constructive field membership (basic)" begin
     R, (x,) = polynomial_ring(Nemo.QQ, ["x"])

     generators = [x^2, x^3]
     to_be_reduced = [x^2, x, 3one(R), zero(R)]

     memberships, remainders, relations_between_tags, tag_to_gen =
         StructuralIdentifiability.check_constructive_field_membership(
             StructuralIdentifiability.RationalFunctionField(generators),
             map(f -> f // one(f), to_be_reduced),
             tag_names = ["T1", "T2"],
         )
     tags = gens(base_ring(parent(first(remainders))))

     @test length(tags) == 2
     @test all(memberships)
     @test map(string, remainders) == ["T1", "T1^2//T2", "3", "0"]
     @test tag_to_gen == Dict(tags[1] => x^2, tags[2] => x^3)
     @test length(relations_between_tags) == 1
     @test string(relations_between_tags[1]) == "T1^3 - T2^2"
end

@testset "RationalFunctionField: constructive field membership" begin
    cases = []

    R, (x, y, z) = Nemo.polynomial_ring(Nemo.QQ, ["x", "y", "z"])

    push!(
        cases,
        Dict(
            :field => RationalFunctionField([[R(1), x + y], [R(1), x * y], [z, (x + y)^2]]),
            :funcs => [
                (x^2 + y^2) // R(1),
                (x^3 + y^3) // (z - x * y),
                R(1) // (z + x + y),
                z // x,
            ],
            :correct => [true, true, true, false],
        ),
    )

    push!(
        cases,
        Dict(
            :field => RationalFunctionField([[
                x + y + z,
                x^2 + y^2 + z^2,
                (x + y + z)^2,
                x^3 + y^3 + z^3,
            ]]),
            :funcs => [x + y + z // 1, x * y * z // 1, x + y + 2 * z // 1, x // (y + z)],
            :correct => [true, true, false, false],
        ),
    )

    push!(
        cases,
        Dict(
            :field => RationalFunctionField([
                x + y + z // 1,
                x * y + y * z + z * x // 1,
                x * y * z // 1,
            ]),
            :funcs => [x^2 + y^2 + z^2, x^6 + y^6 + z^6, x - y + z, x^2 - y^2 + z^2],
            :correct => [true, true, false, false],
        ),
    )

    R, (a, b, c) = QQ["a", "b", "c"]
    F = StructuralIdentifiability.RationalFunctionField([a, b, a + b + c])
    push!(cases, Dict(:field => F, :funcs => [zero(R), one(R)], :correct => [true, true]))
    push!(cases, Dict(:field => F, :funcs => [7a, 9b, 11c], :correct => [true, true, true]))

    F = StructuralIdentifiability.RationalFunctionField([2c, 3b, 5a])
    push!(cases, Dict(:field => F, :funcs => [7a, 9b, 11c], :correct => [true, true, true]))

    s1, s2 = a + b + c, a^2 + b^2 + c^2
    F = StructuralIdentifiability.RationalFunctionField([s1, s2])
    push!(
        cases,
        Dict(
            :field => F,
            :funcs => [a, b + c, a * b + b * c],
            :correct => [false, false, false],
            :relations => [],
        ),
    )
    push!(
        cases,
        Dict(
            :field => F,
            :funcs => [a * b + b * c + a * c, (s1)^8 - (s2)^9 + 89],
            :correct => [true, true],
        ),
    )

    # Takes too much memory and time
    # Example in Section 5 from
    # https://mediatum.ub.tum.de/doc/685465/685465.pdf
    # R, (x1, x2) = QQ["x1", "x2"]
    # g1 = (x1^3 + x1 * x2 - 2) // (x1^2 - x2 - 1)
    # g2 = (x1^2 + x1^2 * x2 + 7) // (x1 - x1^2 * x2^2)
    # g3 = x1^2 + 3x1 * x2
    # g4 = x1 * x2^2 + 5x1 * x2
    # g5 = x1^3 * x2 - x2
    # F = StructuralIdentifiability.RationalFunctionField([g1, g2, g3, g4, g5])
    # push!(cases, Dict(:field => F, :funcs => [x1, x2], :correct => [true, true]))

    R, (a, b, c) = polynomial_ring(QQ, ["a", "b", "c"])
    F = RationalFunctionField([
        (a^2 + b^2) // one(R),
        (a^3 + b^3) // one(R),
        (a^4 + b^4) // one(R),
    ])
    push!(
        cases,
        Dict(
            :field => F,
            :funcs => [
                a // one(R),
                (a + b) // one(R),
                c // a,
                b // a,
                a * b // (a + b),
                c^2 // c + 1,
            ],
            :correct => [false, true, false, false, true, false],
        ),
    )

    R, (T1,) = polynomial_ring(Nemo.QQ, ["T1"])
    push!(
        cases,
        Dict(
            :field => RationalFunctionField([T1^2]),
            :funcs => [T1, T1^2],
            :correct => [false, true]
        )
    )

    R, (T1, t, _t) = polynomial_ring(Nemo.QQ, ["T1", "t", "_t"])
    push!(
        cases,
        Dict(
            :field => RationalFunctionField([T1, t, _t]),
            :funcs => [_t, t, T1 * t * _t],
            :correct => [true, true, true],
        )
    )

    R, (x,) = polynomial_ring(Nemo.QQ, ["x"])
    push!(
        cases,
        Dict(
            :field => RationalFunctionField([(x - 1), R(1) // (x^5 - 1), x]),
            :funcs => [x^4 + x^3 + x^2 + x + 1, x, R(33) // x^2],
            :correct => [true, true, true],
        )
    )
    push!(
        cases,
        Dict(
            :field => RationalFunctionField([(x^10 + x^9 + x^2 + 1) // (x^7 - x^6 - x^3 + 1)]),
            :funcs => [x,],
            :correct => [false],
        )
    )
    push!(
        cases,
        Dict(
            :field => RationalFunctionField([x^2]),
            :funcs => [x, x^42],
            :correct => [false, true],
        )
    )

    R, (x, y, z) = polynomial_ring(Nemo.QQ, ["x", "y", "z"])
    push!(
        cases,
        Dict(
            :field => RationalFunctionField([x^2 + y^2, x^3 + y^3, x^4 + y^4]),
            :funcs => [x * y, x + y, x + y + 1, x + y + z],
            :correct => [true, true, true, false],
        )
    )

    # NOTE: in this case it actually matter to cancel out the gcd after
    # computing the normal forms
    R, (a, b, y, x2, c, x1) = polynomial_ring(Nemo.QQ, ["a", "b", "y", "x2", "c", "x1"])
    push!(
        cases,
        Dict(
            :field => RationalFunctionField([
                x1 // one(R),
                a // one(R),
                (a * c + c^2) // one(R),
                c // x2,
                x2 // (a + b),
            ]),
            :funcs => [
                (a * c + c^2 + x1) // (a * c + c^2),
                (a * c + c^2 + x1) // (a^2 + a * b + a * c + b * c),
                (a * x2 + a * x1 + b * x1) // x2,
            ],
            :correct => [true, true, true],
        )
    )



    for c in cases
        R = poly_ring(c[:field])
        containment, expressions, relations, tag_to_gen = check_constructive_field_membership(c[:field], c[:funcs] .// one(R))
        @test containment == c[:correct]
        for (i, expr) in enumerate(expressions)
            if containment[i]
                @test c[:funcs][i] == StructuralIdentifiability.eval_at_dict(expr, tag_to_gen)
            end
        end
        if :relations in values(c)
            @test Set(relations) == Set([parent_ring_change(p, parent(first(relations))) for p in c[:relations]])
        end
    end
end
