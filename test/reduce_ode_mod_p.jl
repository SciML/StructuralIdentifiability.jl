@testset "Reducing ODE mod p" begin
    ode = @ODEmodel(
        x'(t) = 1 // 2 * x(t) - 5 * y(t),
        y'(t) = (3 // 5 * x(t) + 17 * y(t)) // (y(t) - 1),
        z(t) = y(t)
    )

    ode_red = reduce_ode_mod_p(ode, 17)

    x = str_to_var("x", ode_red.poly_ring)
    y = str_to_var("y", ode_red.poly_ring)
    @test ode_red.x_equations[x] == 9 * x + 12 * y
    @test ode_red.x_equations[y] == 4 * x // (y + 16)
end
