import StructuralIdentifiability: parent_ring_change

@testset "Identifiable functions of parameters" begin
    # For each of the ODEs, 
    test_cases = []

    # Toy models
    # TODO: uncomment when globabl identifiability can handle states
    # ode = @ODEmodel(x'(t) = x(t), y(t) = x(t))
    # T = typeof(x)
    # ident_funcs = Vector{AbstractAlgebra.Generic.Frac{T}}()
    # push!(test_cases, (ode = ode, ident_funcs = ident_funcs))
    #
    ode = @ODEmodel(x'(t) = a * x(t) + u(t), y(t) = x(t))
    ident_funcs = [a // one(a)]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Example 2 from 
    # "On Global Identifiability for Arbitrary Model Parametrizations",
    # DOI: 10.1016/0005-1098(94)90029-9
    ode = @ODEmodel(x1'(t) = Θ * x2(t)^2, x2'(t) = u, y(t) = x1(t))
    ident_funcs = [u^2 * Θ // one(Θ)]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Example 4 from 
    # "On Global Identifiability for Arbitrary Model Parametrizations",
    # DOI: 10.1016/0005-1098(94)90029-9
    ode = @ODEmodel(x'(t) = (-V_m * x(t)) / (k_m + x(t)) + k01 * x(t), y(t) = c * x(t))
    ident_funcs = [k01 // one(k01), k_m // V_m, V_m * c // one(V_m)]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # TODO: uncomment when identifiability can handle zero states 
    # ode = @ODEmodel(
    #     y(t) = a*u(t)
    # )
    # ident_funcs = [(a01 * a12) // one(a01), (a01 + a12 + a21) // one(a01)]
    # push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # 2-compartiment model
    ode = @ODEmodel(
        x0'(t) = -(a01 + a21) * x0(t) + a12 * x1(t),
        x1'(t) = a21 * x0(t) - a12 * x1(t),
        y(t) = x0(t)
    )
    ident_funcs = [(a01 * a12) // one(a01), (a01 + a12 + a21) // one(a01)]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    p = 0.99
    for case in test_cases
        ode = case.ode
        true_ident_funcs = case.ident_funcs
        result_funcs = find_identifiable_functions(ode)
        result_dennums = StructuralIdentifiability.fractions_to_dennums(result_funcs)
        R = parent(result_dennums[1][1])
        println(R)
        println(result_funcs)
        println(true_ident_funcs)
        println(parent(true_ident_funcs[1]))
        true_ident_funcs = map(
            f ->
                parent_ring_change(numerator(f), R) //
                parent_ring_change(denominator(f), R),
            true_ident_funcs,
        )
        true_dennums = StructuralIdentifiability.fractions_to_dennums(true_ident_funcs)
        # Check inclusion <true funcs> in <result funcs>
        inclusion = StructuralIdentifiability.check_field_membership(
            result_dennums,
            true_ident_funcs,
            p,
        )
        @test all(inclusion)
        # Check inclusion <new funcs> in <true funcs>
        inclusion =
            StructuralIdentifiability.check_field_membership(true_dennums, result_funcs, p)
        @test all(inclusion)
    end
end
