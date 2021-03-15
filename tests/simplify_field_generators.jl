@testset "Check identifiability from given IO-equations" begin
    R, (x, y, z) = PolynomialRing(Nemo.QQ, ["x", "y", "z"])

    println(simplify_field_generators([[R(1), x + y, x^2 + y^2]]))

    println(simplify_field_generators([[R(1), (x + y)^2, x^3 + y^3, x^2 + y^2]]))
end
