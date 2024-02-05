cases = []

R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
push!(
    cases,
    Dict(
        :F => RationalFunctionField([[one(R), x + y, x * y]]),
        :funcs => [x // one(R), z // one(R), x^3 - y^3 // one(R), x + z // one(R)],
        :correct => [true, false, true, false],
    ),
)

push!(
    cases,
    Dict(
        :F => RationalFunctionField([[x, y], [y, z]]),
        :funcs => [x // z, (x + y) // z, x // one(R), y // one(R), z // one(R)],
        :correct => [true, true, false, false, false],
    ),
)

@testset "Algebraicity over a field" begin
    for case in cases
        @test check_algebraicity(case[:F], case[:funcs], 0.99) == case[:correct]
    end
end
