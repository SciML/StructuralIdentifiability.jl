@testset "IO-equation of predator-prey model" begin
    # 2-compartiment model
    var_names = [
        "x_0", "x_1", "u",
        "a_01", "a_21", "a_12"
    ]
    R, (x_0, x_1, u, a_01, a_21, a_12) = PolynomialRing(QQ, var_names)
    f = [-(a_01 + a_21) * x_0 + a_12 * x_1 + u, a_21 * x_0 - a_12 * x_1]
    g = x_0

    ode = ODE([x_0, x_1], [a_01, a_21, a_12], f, g)
    find_ioequation(ode, true)

    (x_0, x_1, u, a_01, a_21, a_12, u_1, u_2, x_0_dot, x_1_dot, y_0, y_1, y_2) = gens(ode.poly_ring)
    ind_multi = ode.io_equation // (y_2 + y_1*a_01 + y_1*a_21 + y_1*a_12 + y_0*a_01*a_12 - u*a_12 - u_1)
    @test numerator(ind_multi) == to_base_ring(numerator(ind_multi))
end