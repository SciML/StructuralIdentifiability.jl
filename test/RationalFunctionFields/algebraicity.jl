@testset "Transcendence basis computations and algebraicity checks" begin
    R, (a, b, c) = QQ["a", "b", "c"]

    F = RationalFunctionField([
        (a^2 + b^2) // (a^3 + b^3),
        a * b // (a + b),
        (a + b) // (a^2 + b^3),
    ])
    update_trbasis_info!(F, 0.9)
    @test F.trbasis == [(a^2 + b^2) // (a^3 + b^3), a * b // (a + b)]
    @test F.trbasis_over == [c]
    @test F.trbasis_probability == 0.9

    @test check_algebraicity(F, [a, b, c, a^2 + c^2], 0.99) == [true, true, false, false]

    F = RationalFunctionField([a^10 // a, 1 // c^5])
    @test check_algebraicity(F, [a, b, c], 0.99) == [true, false, true]
    @test F.trbasis == [a^10 // a, 1 // c^5]
    @test F.trbasis_over == [b]
end
