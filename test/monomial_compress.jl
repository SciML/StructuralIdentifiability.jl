@testset "Monomial compression test" begin
    
    R, (v1, v2, v3, v4) = Nemo.PolynomialRing(Nemo.QQ, ["v1", "v2", "v3", "v4"])
    tests = [
        v1 + v2 - 2,
        v4 + v4 * v1 - 3 * v2 * v4,
        (v1 + v2 + v2^2) * (v3 + v4),
        (v1 + v2^2) * (v3^3 + v4) - 7 * (v1 - 3 - v2) * (v3 - v4^2 - v3^2)
    ]

    for t in tests
        a, b = monomial_compress(t, [v1, v2])
        s = 0
        for (x, y) in zip(a, b)
            s += parent_ring_change(x, parent(t)) * parent_ring_change(y, parent(t))
        end
        @test s == t
    end

end
