@testset "Transcendence basis computations and algebraicity checks" begin
    R, (a, b, c) = QQ["a", "b", "c"]

    F = RationalFunctionField([
        (a^3 + b^3) // (a^2 + b^2),
        a * b // (a + b),
        (a + b) // (a^2 + b^3),
    ])
    update_trbasis_info!(F, 0.999)
    @test F.trbasis == [(a^3 + b^3) // (a^2 + b^2), a * b // (a + b)]
    @test F.trbasis_over == [c]
    @test F.trbasis_probability == 0.999

    @test check_algebraicity(F, [a, b, c, a^2 + c^2], 0.99) == [true, true, false, false]
    @test check_algebraicity_modp(F, [a, b, c, a^2 + c^2]) == [true, true, false, false]

    F = RationalFunctionField([a^10 // c, c^5 // one(R)])
    @test check_algebraicity(F, [a, b, c], 0.99) == [true, false, true]
    @test check_algebraicity_modp(F, [a, b, c]) == [true, false, true]
    @test F.trbasis == [a^10 // c, c^5 // one(R)]
    @test F.trbasis_over == [b]

    # Coming from this linear compartment model
    #  linear_compartment_model(
    #  [[2], [3], [4], [5], [1]], # graph
    #  inputs = [1],
    #  outputs = [5],
    #  leaks = [2, 3]
    #  )
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
    # the correct answer verified with SIAN
    correct = [true, false, false, true, true, false, false]
    @test check_algebraicity(F, [x // one(R) for x in gens(R)], 0.999) == correct
    @test check_algebraicity_modp(F, [x // one(R) for x in gens(R)], 2^30 + 3) == correct

    R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
    F = RationalFunctionField([[one(R), x + y, x * y]])
    fs = [x // one(R), z // one(R), x^3 - y^3 // one(R), x + z // one(R)]
    @test check_algebraicity(F, fs, 0.99) == [true, false, true, false]

    F = RationalFunctionField([[x, y], [y, z]])
    fs = [x // z, (x + y) // z, x // one(R), y // one(R), z // one(R)]
    @test check_algebraicity(F, fs, 0.95) == [true, true, false, false, false]
end
