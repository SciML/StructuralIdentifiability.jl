@testset "Identifiable functions of parameters" begin
    # For each of the ODEs, 
    test_cases = []

    # Toy models
    ode = @ODEmodel(x'(t) = x(t), y(t) = x(t))
    T = typeof(x)
    ident_funcs = Vector{AbstractAlgebra.Generic.Frac{T}}()
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))
    #
    ode = @ODEmodel(x'(t) = a * x(t) + u(t), y(t) = x(t))
    ident_funcs = [a // one(a)]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # 2-compartiment model
    ode = @ODEmodel(
        x0'(t) = -(a01 + a21) * x0(t) + a12 * x1(t),
        x1'(t) = a21 * x0(t) - a12 * x1(t),
        y(t) = x0(t)
    )
    ident_funcs = [a01 * a12, a01 + a12 + a21, (a01 + a12 + a21) // (a01 * a12), x0]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    p = 0.99
    for case in test_cases
        ode = case.ode
        true_ident_funcs = case.ident_funcs
        result = find_identifiable_functions(ode)
        # Check inclusion <true funcs> in <new funcs>
        inclusion = StructuralIdentifiability.check_field_membership(result, true_ident_funcs, p)
        @test all(inclusion)
    end
end
