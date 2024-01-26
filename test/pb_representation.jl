@testset "PB-representation - creation" begin
    ode = @ODEmodel(x1'(t) = x2(t), x2'(t) = a * x1(t), y(t) = x1(t))
    ioeqs = find_ioequations(ode)
    pbr = PBRepresentation(ode, ioeqs)
    @test pbr.profile == Dict("y(t)" => 2)
    @test pbr.u_names == Array{String, 1}()
    @test pbr.y_names == ["y(t)"]
    @test pbr.param_names == ["a"]

    ode = @ODEmodel(
        x1'(t) = x1(t),
        x2'(t) = a * x2(t),
        y(t) = x1(t),
        y2(t) = x2(t) + u(t)
    )
    ioeqs = find_ioequations(ode)
    pbr = PBRepresentation(ode, ioeqs)
    @test pbr.profile == Dict("y(t)" => 1, "y2(t)" => 1)
    @test pbr.u_names == ["u(t)"]
    @test pbr.y_names == ["y(t)", "y2(t)"]
    @test pbr.param_names == ["a"]
end
