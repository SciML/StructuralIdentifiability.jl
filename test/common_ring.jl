if GROUP == "All" || GROUP == "Core"
    @testset "Computing common ring for the PB-reduction" begin
        ode = @ODEmodel(x1'(t) = x2(t), x2'(t) = a * x1(t), y(t) = x1(t))
        ioeqs = find_ioequations(ode)
        pbr = PBRepresentation(ode, ioeqs)
        R, (y_2, y_5, c) = Nemo.polynomial_ring(Nemo.QQ, ["y(t)_2", "y(t)_5", "c"])
        p = y_2^2 + c * y_5
        (r, der) = common_ring(p, pbr)
        @test Set(map(var_to_str, gens(r))) ==
              Set(["y(t)_0", "y(t)_1", "y(t)_2", "y(t)_3", "y(t)_4", "y(t)_5", "c", "a"])

        ode = @ODEmodel(
            x1'(t) = x3(t),
            x2'(t) = a * x2(t),
            x3'(t) = x1(t),
            y1(t) = x1(t),
            y2(t) = x2(t) + u(t)
        )
        ioeqs = find_ioequations(ode)
        pbr = PBRepresentation(ode, ioeqs)
        R, (y1_0, y2_3, u_3) =
            Nemo.polynomial_ring(Nemo.QQ, ["y1(t)_0", "y2(t)_3", "u(t)_3"])
        p = y1_0 + y2_3 + u_3
        (r, der) = common_ring(p, pbr)
        @test Set([var_to_str(v) for v in gens(r)]) == Set([
            "y1(t)_0",
            "y1(t)_1",
            "y1(t)_2",
            "y1(t)_3",
            "y1(t)_4",
            "y2(t)_0",
            "y2(t)_1",
            "y2(t)_2",
            "y2(t)_3",
            "u(t)_0",
            "u(t)_1",
            "u(t)_2",
            "u(t)_3",
            "a",
        ])
    end
end
