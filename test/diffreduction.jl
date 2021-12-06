@testset "Differential reduction" begin
    ode = @ODEmodel(
        x'(t) = a * x(t),
        y(t) = x(t)
    )
    pbr = PBRepresentation(ode, find_ioequations(ode))
    R, (y_5, y_4) = Nemo.PolynomialRing(Nemo.QQ, ["y_5", "y_4"])
    res = diffreduce(y_5, pbr)
    @test res == str_to_var("a", parent(res))^5 * str_to_var("y_0", parent(res))
    res = diffreduce(y_4, pbr)
    @test res == str_to_var("a", parent(res))^4 * str_to_var("y_0", parent(res))
    res = diffreduce(y_4^2, pbr)
    @test res == str_to_var("a", parent(res))^8 * str_to_var("y_0", parent(res))^2

    ode = @ODEmodel(
        x1'(t) = x2(t),
        x2'(t) = -x1(t),
        y(t) = x1(t) + u(t)
    )
    pbr = PBRepresentation(ode, find_ioequations(ode))
    R, (y_5, y_4, u_10) = Nemo.PolynomialRing(Nemo.QQ, ["y_5", "y_4", "u_10"])
    res = diffreduce(y_4 + u_10, pbr)
    @test res == str_to_var("u_10", parent(res)) + str_to_var("y_0", parent(res)) - str_to_var("u_0", parent(res)) + str_to_var("u_4", parent(res))

    # Mizuka's example
    R, (y_0, y_1, y_2, y_3, u_0, u_1) = Nemo.PolynomialRing(Nemo.QQ, ["y_0", "y_1", "y_2", "y_3", "u_0", "u_1"])
    pbr = PBRepresentation(["y"], ["u"], Array{String, 1}(), Dict("y" => 2), Dict("y" => 27*y_0^9 + 27*y_0^6*y_1^3 - 54*u_0*y_0^6 + 54*y_0^7 + 54*y_0^6*y_1 - 27*y_0^5*y_1^2 - 27*y_0^4*u_1*y_1^2 + 27*y_0^4*y_1^3 + 27*y_0^4*y_1^2*y_2 + 27*u_0^2*y_0^3 - 54*u_0*y_0^4 - 54*u_0*y_0^3*y_1 + 27*y_0^5 + 54*y_0^4*y_1 + 18*y_0^3*u_1*y_1 + 9*y_0^3*y_1^2 - 18*y_0^3*y_1*y_2 + 9*y_0^2*u_1^2*y_1 - 18*y_0^2*u_1*y_1^2 - 18*y_0^2*u_1*y_1*y_2 + 9*y_0^2*y_1^3 + 18*y_0^2*y_1^2*y_2 + 9*y_0^2*y_1*y_2^2 + 4*y_0^3 - 3*y_0*u_1^2 + 6*y_0*u_1*y_1 + 6*y_0*u_1*y_2 - 3*y_0*y_1^2 - 6*y_0*y_1*y_2 - 3*y_0*y_2^2 - u_1^3 + 3*u_1^2*y_1 + 3*u_1^2*y_2 - 3*u_1*y_1^2 - 6*u_1*y_1*y_2 - 3*u_1*y_2^2 + y_1^3 + 3*y_1^2*y_2 + 3*y_1*y_2^2 + y_2^3))
    p = 100*y_0^4 + 10*y_0^3*y_1 + 410//3*y_0*y_2 + 10*y_0*y_3 - 110*y_1^2 - 10*y_1*y_2
    @time diffreduce(p, pbr)
end
