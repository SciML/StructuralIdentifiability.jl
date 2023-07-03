import StructuralIdentifiability: parent_ring_change

@testset "Identifiable functions of parameters" begin
    # For each ODE system we check the equality (in terms of fields of rational
    # functions) of the true set of identifiable functions and the obtained
    # simplified set
    test_cases = []

    # TODO: uncomment when global identifiability can handle states.
    # ode = @ODEmodel(x'(t) = x(t), y(t) = x(t))
    # T = typeof(x)
    # ident_funcs = Vector{AbstractAlgebra.Generic.Frac{T}}()
    # push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    ode = @ODEmodel(x'(t) = a * x(t) + u(t), y(t) = x(t))
    ident_funcs = [a]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Parameter a is not identifiable, and neither are any combinations thereof.
    ode = @ODEmodel(x1'(t) = x2(t) - a, x2'(t) = x1(t) + a, y(t) = x1(t) + x2(t))
    ident_funcs = Vector{AbstractAlgebra.Generic.Frac{typeof(x1)}}()
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Example 2 from 
    # "On Global Identifiability for Arbitrary Model Parametrizations",
    # DOI: 10.1016/0005-1098(94)90029-9
    ode = @ODEmodel(x1'(t) = Θ * x2(t)^2, x2'(t) = u, y(t) = x1(t))
    # TODO: do we want u^2 Θ or Θ in the output?
    ident_funcs = [u^2 * Θ]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Example 4 from 
    # "On Global Identifiability for Arbitrary Model Parametrizations",
    # DOI: 10.1016/0005-1098(94)90029-9
    ode = @ODEmodel(x'(t) = (-V_m * x(t)) / (k_m + x(t)) + k01 * x(t), y(t) = c * x(t))
    ident_funcs = [k01, k_m // V_m, V_m * c // one(V_m)]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Parameters b and c enter the io equations only as the product b * c.
    ode = @ODEmodel(x'(t) = a * x(t) + b * u(t), y(t) = c * x(t))
    ident_funcs = [b * c // one(b), a // one(a)]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # 2-compartiment model
    ode = @ODEmodel(
        x0'(t) = -(a01 + a21) * x0(t) + a12 * x1(t),
        x1'(t) = a21 * x0(t) - a12 * x1(t),
        y(t) = x0(t)
    )
    ident_funcs = [(a01 * a12), (a01 + a12 + a21)]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # TODO: uncomment when identifiability can handle zero states 
    # ode = @ODEmodel(
    #     y(t) = a*u(t)
    # )
    # ident_funcs = [(a01 * a12), (a01 + a12 + a21)]
    # push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Example 2 from
    # "On Structural Identifiability",
    # DOI: https://doi.org/10.1016/0025-5564(70)90132-X
    #
    # More or less the same 2-compartmental model as the one given above
    ode = @ODEmodel(
        x1'(t) = -(k1 + k2) * x1(t) + k3 * x2(t) + u(t),
        x2'(t) = k2 * x1(t) - (k3 + k4) * x2(t),
        y(t) = x1(t)
    )
    ident_funcs = [(k1 + k2), (k1 + k2 + k3 + k4), ((k1 + k2) * (k3 + k4) - k2 * k3)]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Diagonal with simple spectrum and observable states
    ode = @ODEmodel(
        x1'(t) = λ1 * x1(t) + β1 * u1(t),
        x2'(t) = λ2 * x2(t) + β2 * u2(t),
        x3'(t) = λ3 * x3(t) + β3 * u3(t),
        y(t) = x1(t) + x2(t) + x3(t)
    )
    ident_funcs = [λ1, λ2, λ3, β1, β2, β3]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # 3 compartments:
    #   x1 <--> x2 <--> x3
    # If we observe x1 and control x1, then all parameters are identifiable
    ode = @ODEmodel(
        x1'(t) = -a1 * x1(t) + b1 * x2(t) + u(t),
        x2'(t) = -(a2 + b1) * x2(t) + a1 * x1(t) + b2 * x3(t),
        x3'(t) = -b2 * x3(t) + a2 * x2(t),
        y(t) = x1(t)
    )
    ident_funcs = [a1, a2, b1, b2]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Example 3 from
    # "On Structural Identifiability",
    # DOI: https://doi.org/10.1016/0025-5564(70)90132-X
    # 
    # 3 compartments:
    #   x1 <--> x2 <--> x3
    # If we observe x1 and control x3, then only some functions of parameters
    # are identifiable
    ode = @ODEmodel(
        x1'(t) = -a1 * x1(t) + b1 * x2(t),
        x2'(t) = -(a2 + b1) * x2(t) + a1 * x1(t) + b2 * x3(t),
        x3'(t) = -b2 * x3(t) + a2 * x2(t) + u(t),
        y(t) = x1(t)
    )
    ident_funcs = [b1 * b2, a1 + a2 + b1 + b2, a1 * a2 + a1 * b2 + b1 * b2]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # 
    ode = @ODEmodel(
        x1'(t) = -a1 * x1(t) + b1 * x2(t),
        x2'(t) = -(a2 + b1) * x2(t) + a1 * x1(t) + b2 * x3(t),
        x3'(t) = -b2 * x3(t) + a2 * x2(t) + u(t),
        y(t) = x1(t)
    )
    ident_funcs = [b1 * b2, a1 + a2 + b1 + b2, a1 * a2 + a1 * b2 + b1 * b2]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    p = 0.99
    for case in test_cases
        ode = case.ode
        true_ident_funcs = case.ident_funcs
        result_funcs = find_identifiable_functions(ode)
        result_dennums = StructuralIdentifiability.fractions_to_dennums(result_funcs)
        R = parent(result_dennums[1][1])
        true_ident_funcs = map(f -> f // one(f), true_ident_funcs)
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
