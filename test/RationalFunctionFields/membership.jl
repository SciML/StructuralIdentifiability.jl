@testset "RationalFunctionField: membership" begin
    cases = []

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

    # Example in Section 5 from
    # https://mediatum.ub.tum.de/doc/685465/685465.pdf
    R, (x1, x2) = QQ["x1", "x2"]
    g1 = (x1^3 + x1 * x2 - 2) // (x1^2 - x2 - 1)
    g2 = (x1^2 + x1^2 * x2 + 7) // (x1 - x1^2 * x2^2)
    g3 = x1^2 + 3x1 * x2
    g4 = x1 * x2^2 + 5x1 * x2
    g5 = x1^3 * x2 - x2
    F = StructuralIdentifiability.RationalFunctionField([g1, g2, g3, g4, g5])
    push!(cases, Dict(:field => F, :funcs => [x1, x2], :correct => [true, true]))

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

    #  linear_compartment_model(
    # [[2], [3], [4], [5], [1]], # graph
    # [1], # inputs
    # [5], # outputs
    # [2, 3] # leaks
    # )
    R, (a_2_1, a_3_2, a_4_3, a_5_4, a_1_5, a_0_2, a_0_3) =
        polynomial_ring(QQ, ["a_2_1", "a_3_2", "a_4_3", "a_5_4", "a_1_5", "a_0_2", "a_0_3"])

    F = RationalFunctionField([
        (
            -a_2_1 * a_3_2 * a_4_3 - a_2_1 * a_3_2 * a_5_4 - a_2_1 * a_3_2 * a_1_5 -
            a_2_1 * a_3_2 * a_0_3 - a_2_1 * a_4_3 * a_5_4 - a_2_1 * a_4_3 * a_1_5 -
            a_2_1 * a_4_3 * a_0_2 - a_2_1 * a_5_4 * a_1_5 - a_2_1 * a_5_4 * a_0_2 -
            a_2_1 * a_5_4 * a_0_3 - a_2_1 * a_1_5 * a_0_2 - a_2_1 * a_1_5 * a_0_3 -
            a_2_1 * a_0_2 * a_0_3 - a_3_2 * a_4_3 * a_5_4 - a_3_2 * a_4_3 * a_1_5 -
            a_3_2 * a_5_4 * a_1_5 - a_3_2 * a_5_4 * a_0_3 - a_3_2 * a_1_5 * a_0_3 -
            a_4_3 * a_5_4 * a_1_5 - a_4_3 * a_5_4 * a_0_2 - a_4_3 * a_1_5 * a_0_2 -
            a_5_4 * a_1_5 * a_0_2 - a_5_4 * a_1_5 * a_0_3 - a_5_4 * a_0_2 * a_0_3 -
            a_1_5 * a_0_2 * a_0_3
        ) // (a_2_1 * a_3_2 * a_4_3 * a_5_4),
        (-a_3_2 * a_1_5 * a_0_3 - a_4_3 * a_1_5 * a_0_2 - a_1_5 * a_0_2 * a_0_3) //
        (a_3_2 * a_4_3),
        (-a_2_1 - a_3_2 - a_4_3 - a_5_4 - a_1_5 - a_0_2 - a_0_3) //
        (a_2_1 * a_3_2 * a_4_3 * a_5_4),
        (
            -a_2_1 * a_3_2 * a_4_3 * a_5_4 - a_2_1 * a_3_2 * a_4_3 * a_1_5 -
            a_2_1 * a_3_2 * a_5_4 * a_1_5 - a_2_1 * a_3_2 * a_5_4 * a_0_3 -
            a_2_1 * a_3_2 * a_1_5 * a_0_3 - a_2_1 * a_4_3 * a_5_4 * a_1_5 -
            a_2_1 * a_4_3 * a_5_4 * a_0_2 - a_2_1 * a_4_3 * a_1_5 * a_0_2 -
            a_2_1 * a_5_4 * a_1_5 * a_0_2 - a_2_1 * a_5_4 * a_1_5 * a_0_3 -
            a_2_1 * a_5_4 * a_0_2 * a_0_3 - a_2_1 * a_1_5 * a_0_2 * a_0_3 -
            a_3_2 * a_4_3 * a_5_4 * a_1_5 - a_3_2 * a_5_4 * a_1_5 * a_0_3 -
            a_4_3 * a_5_4 * a_1_5 * a_0_2 - a_5_4 * a_1_5 * a_0_2 * a_0_3
        ) // (a_2_1 * a_3_2 * a_4_3 * a_5_4),
        (
            -a_2_1 * a_3_2 - a_2_1 * a_4_3 - a_2_1 * a_5_4 - a_2_1 * a_1_5 - a_2_1 * a_0_2 - a_2_1 * a_0_3 - a_3_2 * a_4_3 - a_3_2 * a_5_4 -
            a_3_2 * a_1_5 - a_3_2 * a_0_3 - a_4_3 * a_5_4 - a_4_3 * a_1_5 -
            a_4_3 * a_0_2 - a_5_4 * a_1_5 - a_5_4 * a_0_2 - a_5_4 * a_0_3 -
            a_1_5 * a_0_2 - a_1_5 * a_0_3 - a_0_2 * a_0_3
        ) // (a_2_1 * a_3_2 * a_4_3 * a_5_4),
        -1 // (a_2_1 * a_3_2 * a_4_3 * a_5_4),
    ])

    push!(cases, Dict(:field => F, :funcs => gens(R), :correct => [false for _ in gens(R)]))

    # Drug resistance model
    R, (b, br, f, fr, g, ga, ke, mu, s) =
        polynomial_ring(QQ, ["b", "br", "f", "fr", "g", "ga", "ke", "mu", "s"])
    F = RationalFunctionField([
        fr // one(R),
        g + ga + mu // one(R),
        br + g + ga // one(R),
        -br - 2 * g - 2 * ga - mu // one(R),
        ga * s // one(R),
        (-b + br) // br,
        (br + g + ga + mu) // s,
        (-br - g - ga - mu) // s,
        (g * ke + ga * ke + ke * mu) // br,
        (b * g + b * ga + b * mu) // br,
        (-9 * br * ga - 9 * g * ga - 9 * ga^2 - 9 * ga * mu) // (g + ga + mu),
        (-9 * g * s - 9 * ga * s - 9 * mu * s) // (br + g + ga + mu),
        (
            b * br + b * g + 2 * b * ga + b * mu - br^2 - br * f - br * g - 2 * br * ga - br * mu
        ) // br,
        (-br * ke - g * ke - ga * ke - ke * mu) // (br * s),
        (-b * br - b * g - b * ga - b * mu) // (br * s),
        (-b * br - b * g - b * ga - b * mu + br^2 + br * g + br * ga + br * mu) // (br * s),
    ])
    push!(
        cases,
        Dict(
            :field => F,
            :funcs => gens(R),
            :correct => [false, false, false, true, false, false, false, false, false],
        ),
    )

    for c in cases
        @test field_contains(c[:field], c[:funcs], 0.99) == c[:correct]
    end
end
