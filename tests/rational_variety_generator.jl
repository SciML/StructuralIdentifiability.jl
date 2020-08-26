@testset "Rational variety generic point generator" begin
    
    R, (x, y, t) = PolynomialRing(QQ, ["x", "y", "t"])
    
    circle = x^2 + y^2 - 1
    xparam = x * (t^2 + 1) - (1 - t^2)
    yparam = y * (t^2 + 1) - 2 * t
    
    gpg = RationalVarietyPointGenerator([xparam, yparam], [t])
    
    for p in Iterators.take(gpg, 20)
        @test eval_at_dict(circle, p) == 0
    end
    
end
