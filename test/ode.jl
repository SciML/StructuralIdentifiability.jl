@testset "ODE struct" begin
    # adding outputs
    ode = StructuralIdentifiability.@ODEmodel(x'(t) = x(t) + a)
    ode = StructuralIdentifiability.add_outputs(ode, Dict("y" => a * x))
    @test StructuralIdentifiability.AbstractAlgebra.nvars(parent(ode)) == 3

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

@testset "ODE unicode" begin
    ode = StructuralIdentifiability.@ODEmodel(
        🐁'(t) = a * 🐁 - b * 🐁 * 🦉,
        🦉'(t) = c * 🦉 + d * 🐁 * 🦉,
        y(t) = 🐁
    )
    println(ode)
    res = StructuralIdentifiability.assess_identifiability(ode)
    println(res)
    @test res == Dict(
        a => :globally,
        b => :nonidentifiable,
        c => :globally,
        d => :globally,
        🐁 => :globally,
        🦉 => :nonidentifiable,
    )

    ode = StructuralIdentifiability.@ODEmodel(
        ⬜'(t) = a⬜ * ⬜ * 🐁b🦉c,
        🐁b🦉c'(t) = 🐁b🦉c,
        🐁y🐁(t) = ⬜
    )
    println(ode)
    StructuralIdentifiability.assess_identifiability(ode)
    StructuralIdentifiability.find_identifiable_functions(ode)
    StructuralIdentifiability.reparametrize_global(ode)
end
