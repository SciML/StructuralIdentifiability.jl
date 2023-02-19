@testset "Set parameter values" begin
    sys = @ODEmodel(x'(t) = a * x(t) + a^2, y(t) = x(t))

    new_sys = set_parameter_values(sys, Dict(a => Nemo.QQ(1, 2)))
    @test evaluate(first(values(new_sys.x_equations)), [2, 2]) == Nemo.QQ(5, 4)

    new_sys = set_parameter_values(sys, Dict(a => Nemo.QQ(1, 20)))
    @test evaluate(first(values(new_sys.x_equations)), [5, 5]) == Nemo.QQ(101, 400)

    new_sys = set_parameter_values(sys, Dict(a => 0.1))
    @test evaluate(first(values(new_sys.x_equations)), [2, 2]) == Nemo.QQ(21, 100)

    new_sys = set_parameter_values(sys, Dict(a => 2.1))
    @test evaluate(first(values(new_sys.x_equations)), [2, 2]) == Nemo.QQ(861, 100)

    sys = @ODEmodel(x'(t) = a + b + c, y(t) = x(t))
    new_sys = set_parameter_values(sys, Dict(a => 2.1, b => 5, c => 3 // 10))
    @test evaluate(first(values(new_sys.x_equations)), [2, 2]) == Nemo.QQ(74, 10)
end
