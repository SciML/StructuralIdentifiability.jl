@testset "PB-representation - creation" begin
    ode = @ODEmodel(
        x1'(t) = x3(t),
        x2'(t) = a * x2(t),
        x3'(t) = x1(t),
        y1(t) = x1(t),
        y2(t) = x2(t) + u(t)
    )
    ioeqs = find_ioequations(ode)
    println("IOEQS: ", ioeqs)
    pbr = PBRepresentation(ode, ioeqs)

    R, (y1_0, y1_1, y1_2, y2_0, y2_1, y2_2) =
        Nemo.PolynomialRing(Nemo.QQ, ["y1_0", "y1_1", "y1_2", "y2_0", "y2_1", "y2_2"])
    @test find_leader([a, x1, x2, u], pbr) == nothing
    @test find_leader([a, x1, y1_0, x2], pbr) == y1_0
    @test find_leader([y1_0, y1_1, y2_0, y2_1], pbr) == y2_1
    l = find_leader([y1_2, y2_1], pbr)
    @test (pbr.y_names[1] == "y1" && l == y2_1) || (pbr.y_names[1] == "y2" && l == y1_2)
    @test find_leader([y1_2, y2_0], pbr) == y1_2
end
