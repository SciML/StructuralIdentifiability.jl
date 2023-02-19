@testset "Decomposing derivative" begin
    cases = [
        ["yy_11", ["y", "yy", "yy_"], ("yy", 11)],
        ["xx_xx_xx_0", ["xx", "x", "xx_xx_xx"], ("xx_xx_xx", 0)],
        ["abc154f", ["ab", "abc"], nothing],
        ["c_1542673", ["a", "b", "c"], ("c", 1542673)],
        ["a", ["a"], nothing],
    ]
    for c in cases
        @test decompose_derivative(c[1], c[2]) == c[3]
    end
end
