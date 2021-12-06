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
    R, (y_0, y_1, y_2, y_3, u_0, u_1, u_2) = Nemo.PolynomialRing(Nemo.QQ, ["y_0", "y_1", "y_2", "y_3", "u_0", "u_1", "u_2"])
    pbr = PBRepresentation(["y"], ["u"], Array{String, 1}(), Dict("y" => 2), Dict("y" => 27*y_0^9 + 27*y_0^6*y_1^3 - 54*u_0*y_0^6 + 54*y_0^7 + 54*y_0^6*y_1 - 27*y_0^5*y_1^2 - 27*y_0^4*u_1*y_1^2 + 27*y_0^4*y_1^3 + 27*y_0^4*y_1^2*y_2 + 27*u_0^2*y_0^3 - 54*u_0*y_0^4 - 54*u_0*y_0^3*y_1 + 27*y_0^5 + 54*y_0^4*y_1 + 18*y_0^3*u_1*y_1 + 9*y_0^3*y_1^2 - 18*y_0^3*y_1*y_2 + 9*y_0^2*u_1^2*y_1 - 18*y_0^2*u_1*y_1^2 - 18*y_0^2*u_1*y_1*y_2 + 9*y_0^2*y_1^3 + 18*y_0^2*y_1^2*y_2 + 9*y_0^2*y_1*y_2^2 + 4*y_0^3 - 3*y_0*u_1^2 + 6*y_0*u_1*y_1 + 6*y_0*u_1*y_2 - 3*y_0*y_1^2 - 6*y_0*y_1*y_2 - 3*y_0*y_2^2 - u_1^3 + 3*u_1^2*y_1 + 3*u_1^2*y_2 - 3*u_1*y_1^2 - 6*u_1*y_1*y_2 - 3*u_1*y_2^2 + y_1^3 + 3*y_1^2*y_2 + 3*y_1*y_2^2 + y_2^3))
    p = 100*y_0^4 + 10*y_0^3*y_1 + 410//3*y_0*y_2 + 10*y_0*y_3 - 110*y_1^2 - 10*y_1*y_2
    @time res = diffreduce(p, pbr)
    expected = parent_ring_change(-100*y_2^2*y_1^2+70*y_2^2*y_0^4+380//3*y_2^2*y_0^2-200*y_2*y_1^3-180*y_2*y_0^7-380*y_2*y_0^5-270*y_1^4*y_0^6-1080*y_1^4*y_0^4-630*y_1^4*y_0^2+810*y_1^3*y_0^9-2520*y_1^3*y_0^7-2910*y_1^3*y_0^5-440*y_1^3*y_0^3+220//3*y_1^3*y_0+90*y_1^2*y_0^8+3690*y_1^2*y_0^6-1330*y_1^2*y_0^4+380*y_1^2*y_0^2+1080*y_1*y_0^9-6540*y_1*y_0^7-7220*y_1*y_0^5+810*u_0^2*y_0^6-3420*u_0^2*y_0^4-1620*u_0*y_0^9+5220*u_0*y_0^7+6840*u_0*y_0^5-10*u_1^3*y_1-30*u_1^3*y_0^3+380//3*u_1^3*y_0-80*u_1^2*y_1^2+10*u_1^2*y_0^4+380*u_1^2*y_0^2+190*u_1*y_1^3+180*u_1*y_0^7+380*u_1*y_0^5-3300*y_0^6-1520//3*y_0^4-20*u_2*y_2*u_1*y_0+60*u_2*y_2*y_1*y_0^3+20*u_2*y_2*y_1*y_0-60*u_2*u_1*y_1*y_0^3-20*u_2*u_1*y_1*y_0-360*y_2*u_1*y_1*y_0^5+1380*y_2*u_1*y_1*y_0^3+1580//3*y_2*u_1*y_1*y_0-30*y_2^2*y_1^2*y_0^2+90*y_2^2*y_1*y_0^5-340*y_2^2*y_1*y_0^3-380//3*y_2^2*y_1*y_0-180*y_2*y_1^3*y_0^4-660*y_2*y_1^3*y_0^2+540*y_2*y_1^2*y_0^7-1860*y_2*y_1^2*y_0^5-1380*y_2*y_1^2*y_0^3-160//3*y_2*y_1^2*y_0+240*y_2*y_1*y_0^6+1400*y_2*y_1*y_0^4+1520//3*y_2*y_1*y_0^2+10*u_2*y_2^2*y_0-20*u_2*y_2*y_0^2+10*u_2*u_1^2*y_0+20*u_2*u_1*y_0^2+90*u_2*y_1^2*y_0^5+60*u_2*y_1^2*y_0^3+10*u_2*y_1^2*y_0-60*u_2*y_1*y_0^4-20*u_2*y_1*y_0^2-10*y_2^2*u_1*y_1-30*y_2^2*u_1*y_0^3+380//3*y_2^2*u_1*y_0+20*y_2*u_1^2*y_1+60*y_2*u_1^2*y_0^3-760//3*y_2*u_1^2*y_0+180*y_2*u_1*y_1^2-80*y_2*u_1*y_0^4-1520//3*y_2*u_1*y_0^2+180*y_2*u_0*y_0^4+30*u_1^2*y_1^2*y_0^2+270*u_1^2*y_1*y_0^5-1040*u_1^2*y_1*y_0^3-400*u_1^2*y_1*y_0+90*u_1*y_1^3*y_0^4+600*u_1*y_1^3*y_0^2-810*u_1*y_1^2*y_0^7+2820*u_1*y_1^2*y_0^5+2170*u_1*y_1^2*y_0^3+200*u_1*y_1^2*y_0-60*u_1*y_1*y_0^6-2100*u_1*y_1*y_0^4-760*u_1*y_1*y_0^2-180*u_1*u_0*y_0^4-1080*y_1*u_0*y_0^6+7020*y_1*u_0*y_0^4-6030*y_0^8+810*y_0^12-1800*y_0^10-100*y_1^4, parent(res))
    @test total_degree(divexact(res, expected)) == 0
end
