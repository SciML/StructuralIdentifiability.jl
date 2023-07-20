# import StructuralIdentifiability: parent_ring_change
# using AbstractAlgebra

# TODO: verify that Maple returns the same
@testset "Identifiable functions of parameters" begin
    # For each ODE system we check the equality (in terms of fields of rational
    # functions) of the true set of identifiable functions and the obtained
    # simplified set
    test_cases = []

    ode = StructuralIdentifiability.@ODEmodel(x'(t) = a * x(t) + u(t), y(t) = x(t))
    ident_funcs = [a]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Parameter a is not identifiable, and neither are any combinations thereof.
    ode = StructuralIdentifiability.@ODEmodel(
        x1'(t) = x2(t) - a,
        x2'(t) = x1(t) + a,
        y(t) = x1(t) + x2(t)
    )
    ident_funcs =
        Vector{StructuralIdentifiability.AbstractAlgebra.Generic.Frac{typeof(x1)}}()

    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Example 2 from 
    # "On Global Identifiability for Arbitrary Model Parametrizations",
    # DOI: 10.1016/0005-1098(94)90029-9
    ode =
        StructuralIdentifiability.@ODEmodel(x1'(t) = Θ * x2(t)^2, x2'(t) = u, y(t) = x1(t))
    # TODO: do we want u^2 Θ or Θ in the output?
    ident_funcs = [u^2 * Θ]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Example 4 from 
    # "On Global Identifiability for Arbitrary Model Parametrizations",
    # DOI: 10.1016/0005-1098(94)90029-9
    ode = StructuralIdentifiability.@ODEmodel(
        x'(t) = (-V_m * x(t)) / (k_m + x(t)) + k01 * x(t),
        y(t) = c * x(t)
    )
    ident_funcs = [k01, k_m // V_m, V_m * c]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Parameters b and c enter the io equations only as the product b * c.
    ode = StructuralIdentifiability.@ODEmodel(x'(t) = a * x(t) + b * u(t), y(t) = c * x(t))
    ident_funcs = [b * c, a]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # 2-compartiment model
    ode = StructuralIdentifiability.@ODEmodel(
        x0'(t) = -(a01 + a21) * x0(t) + a12 * x1(t),
        x1'(t) = a21 * x0(t) - a12 * x1(t),
        y(t) = x0(t)
    )
    ident_funcs = [(a01 * a12), (a01 + a12 + a21)]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # TODO: uncomment when identifiability can handle nos states 
    # ode = StructuralIdentifiability.@ODEmodel(
    #     y(t) = a*u(t)
    # )
    # ident_funcs = [(a01 * a12), (a01 + a12 + a21)]
    # push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Example 2 from
    # "On Structural Identifiability",
    # DOI: https://doi.org/10.1016/0025-5564(70)90132-X
    #
    # More or less the same 2-compartmental model as the one given above
    ode = StructuralIdentifiability.@ODEmodel(
        x1'(t) = -(k1 + k2) * x1(t) + k3 * x2(t) + u(t),
        x2'(t) = k2 * x1(t) - (k3 + k4) * x2(t),
        y(t) = x1(t)
    )
    ident_funcs = [(k1 + k2), (k1 + k2 + k3 + k4), ((k1 + k2) * (k3 + k4) - k2 * k3)]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Diagonal with simple spectrum and observable states
    ode = StructuralIdentifiability.@ODEmodel(
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
    ode = StructuralIdentifiability.@ODEmodel(
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
    ode = StructuralIdentifiability.@ODEmodel(
        x1'(t) = -a1 * x1(t) + b1 * x2(t),
        x2'(t) = -(a2 + b1) * x2(t) + a1 * x1(t) + b2 * x3(t),
        x3'(t) = -b2 * x3(t) + a2 * x2(t) + u(t),
        y(t) = x1(t)
    )
    ident_funcs = [b1 * b2, a1 + a2 + b1 + b2, a1 * a2 + a1 * b2]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Example 3 from
    # "On the identifiability and distinguishability of nonlinear parametric
    # models",
    # DOI: https://doi.org/10.1016/0378-4754(95)00123-9
    ode = StructuralIdentifiability.@ODEmodel(
        x1'(t) = p1 * x1^2 + p2 * x1 * x2,
        x2'(t) = p3 * x1^2 + p4 * x1 * x2,
        y(t) = x1
    )
    ident_funcs = [p1 + p4, p2 * p3 - p4 * p1]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Goowdin oscillator
    ode = StructuralIdentifiability.@ODEmodel(
        x1'(t) = -b * x1(t) + 1 / (c + x4(t)),
        x2'(t) = alpha * x1(t) - beta * x2(t),
        x3'(t) = gama * x2(t) - delta * x3(t),
        x4'(t) = sigma * x4(t) * (gama * x2(t) - delta * x3(t)) / x3(t),
        y(t) = x1(t)
    )
    ident_funcs = [sigma, delta + beta, c, b, delta * beta]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # SIRS forced
    ode = StructuralIdentifiability.@ODEmodel(
        s'(t) = mu - mu * s(t) - b0 * (1 + b1 * x1(t)) * i(t) * s(t) + g * r(t),
        i'(t) = b0 * (1 + b1 * x1(t)) * i(t) * s(t) - (nu + mu) * i(t),
        r'(t) = nu * i(t) - (mu + g) * r(t),
        x1'(t) = -M * x2(t),
        x2'(t) = M * x1(t),
        y1(t) = i(t),
        y2(t) = r(t)
    )
    ident_funcs = [g, mu, b0, nu, M^2]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # SEIR_1_io
    ode = StructuralIdentifiability.@ODEmodel(
        S'(t) = -beta * S(t) * I(t),
        E'(t) = beta * S(t) * I(t) - v * E(t),
        I'(t) = v * E(t) - psi * I(t) - (1 - psi) * gamma * I(t),
        R'(t) = gamma * Q(t) + (1 - psi) * gamma * I(t),
        Q'(t) = -gamma * Q(t) + psi * I(t),
        y1(t) = Q(t)
    )
    ident_funcs = [gamma, beta // psi, gamma * psi - v - psi, gamma * psi - v * psi]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Bilirubin2_io.
    # Regression test: failed before, as the total degrees were being estimated
    # incorrectly in the interpolation
    ode = StructuralIdentifiability.@ODEmodel(
        x1'(t) =
            -(k21 + k31 + k41 + k01) * x1(t) +
            k12 * x2(t) +
            k13 * x3(t) +
            k14 * x4(t) +
            u(t),
        x2'(t) = k21 * x1(t) - k12 * x2(t),
        x3'(t) = k31 * x1(t) - k13 * x3(t),
        x4'(t) = k41 * x1(t) - k14 * x4(t),
        y1(t) = x1(t)
    )
    # TODO
    # ident_funcs = [k01, k31 + k21 + k41, k31 * k21 * k41, k31 * k21 + k31 * k41 + k21 * k41]
    # push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Biohydrogenation_io
    ode = StructuralIdentifiability.@ODEmodel(
        x5'(t) =
            (
                k5 * k8 * x4(t) + k5 * x6(t) * x4(t) + k5 * x5(t) * x4(t) -
                k6 * x5(t) * k7 - x5(t) * k7 * x4(t)
            ) // (
                k8 * k6 +
                k8 * x4(t) +
                k6 * x6(t) +
                k6 * x5(t) +
                x6(t) * x4(t) +
                x5(t) * x4(t)
            ),
        x7'(t) = (k9 * k10 * x6(t) - k9 * x6(t)^2) // k10,
        x4'(t) = (-k5 * x4(t)) // (k6 + x4(t)),
        x6'(t) =
            (
                -k8 * k9 * k10 * x6(t) + k8 * k9 * x6(t)^2 - k9 * k10 * x6(t)^2 -
                k9 * k10 * x6(t) * x5(t) +
                k9 * x6(t)^3 +
                k9 * x6(t)^2 * x5(t) +
                k10 * x5(t) * k7
            ) // (k8 * k10 + k10 * x6(t) + k10 * x5(t)),
        y1(t) = x4(t),
        y2(t) = x5(t)
    )
    # TODO: simplify?
    #   k9 // k10, k10^2
    # into
    #   k9 * k10, k10^2
    ident_funcs = [k7, k6, k5, k9 // k10, k10^2, k8 + 1 // 2 * k10]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # SLIQR
    ode = StructuralIdentifiability.@ODEmodel(
        S'(t) = -b * In(t) * S(t) * Ninv - S(t) * Ninv * u(t),
        In'(t) = -In(t) * g + s * Q(t) + a * L(t),
        L'(t) = b * In(t) * S(t) * Ninv - a * L(t),
        Q'(t) = -e * In(t) * g + In(t) * g - s * Q(t),
        y(t) = In(t) * Ninv
    )
    ident_funcs = [
        g + a,
        s + g + a,
        s,
        Ninv,
        b,
        (e * s - s + a) // (e * s^2 * g - s^2 * g - s^2 * a + s * g * a + s * a^2),
        e * s * g + s * a + g * a,
        (e * s^2 + e * s * g - s^2 - s * g + g * a + a^2) //
        (e * s^2 * g - s^2 * g - s^2 * a + s * g * a + s * a^2),
        e * s * g * a,
        2 * e * Ninv * s * g + 2 * Ninv * s * a + 2 * Ninv * g * a,
    ]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # St.
    # Regression test:
    # Failed before, as the degrees of Groebner basis were too large
    ode = StructuralIdentifiability.@ODEmodel(
        S'(t) = -e * S(t) - S(t) * d * W(t) + S(t) * r - S(t) * a * W(t) + R(t) * g,
        R'(t) = e * S(t) + rR * R(t) + S(t) * a * W(t) - dr * R(t) * W(t) - R(t) * g,
        W'(t) = T * Dd - W(t) * Dd,
        y1(t) = S(t) + R(t),
        y2(t) = T
    )
    ident_funcs = [
        T,
        Dd,
        e - rR + dr * T + d * T + g - r + a * T,
        (2 * rR * d - 2 * dr * r) // (dr - d),
        (dr^2 + d^2 + 2 * d * a + a^2) // (dr * d + dr * a),
        (e * dr - e * d + rR * a + dr * g - d * g - r * a) // (dr - d),
        (e * dr^2 - e * dr * d + rR * dr * a + dr * d * g - dr * r * a - d^2 * g) //
        (dr^2 + dr * a - d^2 - d * a),
    ]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    p = 0.99
    for case in test_cases
        ode = case.ode
        true_ident_funcs = case.ident_funcs
        result_funcs = StructuralIdentifiability.find_identifiable_functions(ode)
        result_funcs_non_simplified =
            StructuralIdentifiability.find_identifiable_functions(ode, simplify = false)

        if isempty(true_ident_funcs)
            @test isempty(result_funcs)
            continue
        end
        R = parent(numerator(result_funcs[1]))
        true_ident_funcs = map(f -> f // one(f), true_ident_funcs)
        true_ident_funcs = map(
            f ->
                StructuralIdentifiability.parent_ring_change(numerator(f), R) //
                StructuralIdentifiability.parent_ring_change(denominator(f), R),
            true_ident_funcs,
        )

        # Check inclusion <true funcs> in <result funcs>
        inclusion = StructuralIdentifiability.check_field_membership(
            result_funcs,
            true_ident_funcs,
            p,
        )
        @test all(inclusion)

        # Check inclusion <result funcs> in <result funcs non simplified>
        # inclusion = StructuralIdentifiability.check_field_membership(
        #     result_funcs_non_simplified,
        #     result_funcs,
        #     p,
        # )
        # @test all(inclusion)

        # Check inclusion <result funcs non simplified> in <true funcs>
        # inclusion = StructuralIdentifiability.check_field_membership(
        #     true_ident_funcs,
        #     result_funcs_non_simplified,
        #     p,
        # )
        # @test all(inclusion)

        # Check inclusion <result funcs> in <true funcs>
        inclusion = StructuralIdentifiability.check_field_membership(
            true_ident_funcs,
            result_funcs,
            p,
        )
        @test all(inclusion)
    end
end

@testset "Identifiable functions of states" begin
    test_cases = []

    # TODO: uncomment when this is handled!
    # ode = StructuralIdentifiability.@ODEmodel(x'(t) = x(t), y(t) = x(t))
    # T = typeof(x)
    # ident_funcs = [x]
    # push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    ode = StructuralIdentifiability.@ODEmodel(x'(t) = a * x(t) + u(t), y(t) = x(t))
    ident_funcs = [a, x]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    ode = StructuralIdentifiability.@ODEmodel(
        x1'(t) = a * x1(t) - b * x1(t) * x2(t),
        x2'(t) = -c * x2(t) + d * x1(t) * x2(t),
        y(t) = x1(t)
    )
    ident_funcs = [x1, c, d, a, b*x2]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    ode = StructuralIdentifiability.@ODEmodel(
        x1'(t) = a * x1(t) - b * x1(t) * x2(t),
        x2'(t) = -c * x2(t) + d * x1(t) * x2(t),
        y(t) = x1(t)
    )
    ident_funcs = [x1, c, d, a, b*x2]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    # Diagonal with simple spectrum and observable states
    ode = StructuralIdentifiability.@ODEmodel(
        x1'(t) = λ1 * x1(t) + β1 * u1(t),
        x2'(t) = λ2 * x2(t) + β2 * u2(t),
        x3'(t) = λ3 * x3(t) + β3 * u3(t),
        y(t) = x1(t) + x2(t) + x3(t)
    )
    ident_funcs = [λ1, λ2, λ3, β1, β2, β3, x1, x2, x3]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    ode = StructuralIdentifiability.@ODEmodel(
        x1'(t) = x1(t) + Θ * x2(t),
        x2'(t) = 0,
        y(t) = x1(t)
    )
    ident_funcs = [x1, x2 * Θ]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    ode = StructuralIdentifiability.@ODEmodel(
        x1'(t) = α * x2(t),
        x2'(t) = x3(t),
        x3'(t) = C,
        y(t) = x1(t)
    )
    ident_funcs = [x1, α * x2, α * x3, α * C]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    ode = StructuralIdentifiability.@ODEmodel(
        x1'(t) = α*(x1 - x2),
        x2'(t) = α*(x1 + x2),
        y(t) = (x1^2 + x2^2) // 2,
    )
    ident_funcs = [α, x1^2 + x2^2]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    ode = StructuralIdentifiability.@ODEmodel(x'(t) = a * x(t) + b * u(t), y(t) = c * x(t))
    ident_funcs = [b * c, a, x//b]
    push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

    p = 0.99
    for case in test_cases
        ode = case.ode
        true_ident_funcs = case.ident_funcs
        result_funcs =
            StructuralIdentifiability.find_identifiable_functions(ode, with_states = true)

        if isempty(true_ident_funcs)
            @test isempty(result_funcs)
            continue
        end
        R = parent(numerator(result_funcs[1]))
        true_ident_funcs = map(f -> f // one(f), true_ident_funcs)
        true_ident_funcs = map(
            f ->
                StructuralIdentifiability.parent_ring_change(numerator(f), R) //
                StructuralIdentifiability.parent_ring_change(denominator(f), R),
            true_ident_funcs,
        )

        # Check inclusion <true funcs> in <result funcs>
        inclusion = StructuralIdentifiability.check_field_membership(
            result_funcs,
            true_ident_funcs,
            p,
        )
        @test all(inclusion)

        # Check inclusion <result funcs> in <true funcs>
        inclusion = StructuralIdentifiability.check_field_membership(
            true_ident_funcs,
            result_funcs,
            p,
        )
        @test all(inclusion)
    end
end
