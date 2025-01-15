@testset "Rational function comparison" begin
    R, (a, b, c) = QQ["a", "b", "c"]

    @test rational_function_cmp(a, b * c)
    @test rational_function_cmp(a, b + c)
    @test rational_function_cmp(a * b, a // b)
    @test rational_function_cmp(a // b, b // a)
    @test rational_function_cmp(a // b, a * b + b * c + c^2)
    @test rational_function_cmp(a^2, a^3)
    @test rational_function_cmp((a^2 + a * b + c) // (a - b), (a - b) // (a^2 + a * b + c))
end
