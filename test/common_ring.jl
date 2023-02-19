@testset "Computing common ring for the PB-reduction" begin
    ode = @ODEmodel(x1'(t) = x2(t), x2'(t) = a * x1(t), y(t) = x1(t))
    ioeqs = find_ioequations(ode)
    pbr = PBRepresentation(ode, ioeqs)
    R, (y_2, y_5, c) = Nemo.PolynomialRing(Nemo.QQ, ["y_2", "y_5", "c"])
    p = y_2^2 + c * y_5
    (r, der) = common_ring(p, pbr)
    @test Set(map(var_to_str, gens(r))) ==
          Set(["y_0", "y_1", "y_2", "y_3", "y_4", "y_5", "c", "a"])

    ode = @ODEmodel(
        x1'(t) = x3(t),
        x2'(t) = a * x2(t),
        x3'(t) = x1(t),
        y1(t) = x1(t),
        y2(t) = x2(t) + u(t)
    )
    ioeqs = find_ioequations(ode)
    pbr = PBRepresentation(ode, ioeqs)
    R, (y1_0, y2_3, u_3) = Nemo.PolynomialRing(Nemo.QQ, ["y1_0", "y2_3", "u_3"])
    p = y1_0 + y2_3 + u_3
    (r, der) = common_ring(p, pbr)
    @test Set([var_to_str(v) for v in gens(r)]) == Set([
        "y1_0",
        "y1_1",
        "y1_2",
        "y1_3",
        "y1_4",
        "y2_0",
        "y2_1",
        "y2_2",
        "y2_3",
        "u_0",
        "u_1",
        "u_2",
        "u_3",
        "a",
    ])
end
