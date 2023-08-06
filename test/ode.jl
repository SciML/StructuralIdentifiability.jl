import AbstractAlgebra

@testset "ODE struct" begin
    # adding outputs
    ode = StructuralIdentifiability.@ODEmodel(x'(t) = x(t) + a)
    ode = StructuralIdentifiability.add_outputs(ode, Dict("y" => a * x))
    @test AbstractAlgebra.nvars(parent(ode)) == 3

    # two or more equations on the same variable
    @test_throws DomainError StructuralIdentifiability.@ODEmodel(
        x1'(t) = x1(t) + a,
        x1'(t) = x1(t) + b,
        y(t) = x1
    )
    @test_throws DomainError StructuralIdentifiability.@ODEmodel(
        x1'(t) = x1(t),
        y(t) = x1,
        y(t) = x1
    )
end
