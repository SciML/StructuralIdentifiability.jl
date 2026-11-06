include(joinpath(@__DIR__, "..", "shared", "test_setup.jl"))

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

    # output appearing in rhs
    @test_throws ArgumentError StructuralIdentifiability.@ODEmodel(
        x1'(t) = a + y(t),
        y(t) = x1(t)
    )
    @test_throws ArgumentError StructuralIdentifiability.@ODEmodel(
        x1'(t) = a,
        y(t) = x1(t) + y(t)
    )
end

# Some nasty corner cases
@testset "Reparametrization names" begin
    suggest_names = StructuralIdentifiability.suggest_names
    ode = @ODEmodel(
        x10'(t) = a1 * x10(t) + b * z(t) + x̃10(t),
        z'(t) = 0,
        X3(t) = x10(t),
    )
    @test suggest_names(ode, ode.x_vars, ode.parameters) ==
        string.(vcat(ode.x_vars, ode.parameters, ode.u_vars, ode.y_vars))
    
    # Preserve indices and reserve names of parameters, inputs, and outputs.
    @test suggest_names(
        ode, [a1 * x10, b * x10, x10 + z, x10] .// one(parent(ode)),
        [a1 + b, a1] .// one(parent(ode)),
    ) == ["x̃10_2(t)", "x̃10_3(t)", "X3_2(t)", "x10(t)", "a1_2", "a1", "x̃10(t)", "X3(t)"]
    @test suggest_names(ode, [x10], [a1, a1]) ==
        ["x10(t)", "a1", "a1_2", "x̃10(t)", "X3(t)"]

    # An unchanged name always wins.
    ode = @ODEmodel(n'(t) = a * n(t), ñ'(t) = 0, y(t) = n(t))
    for states in ([a * n, ñ], [ñ, a * n])
        names = suggest_names(ode, states, [a])
        @test Dict(names[1:2] .=> states) == Dict("ñ(t)" => ñ, "ñ_2(t)" => a * n)
    end

    R, (x0, a, y) = polynomial_ring(QQ, ["x0", "a", "y"])
    ode = ODE{typeof(x0)}(Dict(x0 => a * x0), Dict(y => x0), typeof(x0)[])
    @test suggest_names(ode, [a * x0, 2 * a * x0], [a]) == ["x̃0", "x̃0_2", "a", "y"]

    ode = @ODEmodel(x'(t) = (k8 + k9 + k10) * x(t), y(t) = x(t))
    @test suggest_names(ode, [x], [k9^2]) == ["x(t)", "k̃9", "y(t)"]
    @test suggest_names(ode, [x], [k9^2, 1 // k9, k9]) == ["x(t)", "k̃9", "k̃9_2", "k9", "y(t)"]
    @test suggest_names(ode, [x], [k10 + 2 * k8, k10 * k9]) == ["x(t)", "a1", "a2", "y(t)"]
    @test suggest_names(ode, empty(ode.x_vars), [k9^2]) == ["k̃9", "y(t)"]
end

@testset "Rename variables" begin
    ode = @ODEmodel(
        x'(t) = (a * x(t) + u(t)) // (b + x(t)),
        z'(t) = x(t) - z(t),
        y(t) = x(t) + z(t),
    )
    original = sprint(show, ode)
    @test sprint(show, rename_variables(ode, Dict{String, String}())) == original
    for transform in (
            Dict("x" => "χ", "a" => "α", "u" => "signal", "y" => "output"),
            Dict("x" => "ignored", "x(t)" => "χ(t)", "a" => "α", "u(t)" => "signal(t)", "y(t)" => "output(t)"),
        )
        renamed = rename_variables(ode, transform)
        χ, z, α, b, signal, output = str_to_var.(
            ["χ(t)", "z(t)", "α", "b", "signal(t)", "output(t)"], Ref(parent(renamed)),
        )
        @test renamed.x_vars == [χ, z]
        @test renamed.u_vars == [signal]
        @test renamed.y_vars == [output]
        @test Set(string.(renamed.parameters)) == Set(["α", "b"])
        @test renamed.x_equations[χ] == (α * χ + signal) // (b + χ)
        @test renamed.x_equations[z] == χ - z
        @test renamed.y_equations[output] == χ + z
        restored = rename_variables(renamed, Dict("χ" => "x", "α" => "a", "signal" => "u", "output" => "y"))
        @test sprint(show, restored) == original
        @test sprint(show, ode) == original
    end
end

@testset "ODE/DDE unicode" begin
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

    dde = StructuralIdentifiability.@DDSmodel(
        🐁(t + 1) = a * 🐁(t) - b * 🐁(t) * 🦉(t),
        🦉(t + 1) = c * 🦉(t) + d * 🐁(t) * 🦉(t),
        y(t) = 🐁(t)
    )
    println(dde)
    StructuralIdentifiability.assess_local_identifiability(dde)

    # baseod on https://github.com/SciML/StructuralIdentifiability.jl/issues/455
    ode = @ODEmodel(x₂'(t) = x₂(t) + a₂, y₂(t) = x₂(t) * u₂(t))
    StructuralIdentifiability.assess_identifiability(ode)
end
