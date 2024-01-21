if GROUP == "All" || GROUP == "Core"
    ≡(A, B) = A ⊆ B && B ⊆ A

    function test_reparametrization(old_ode, new_ode, var_mapping, implicit_relations)
        # NOTE: should this be strengthened to == ?
        # @test map(string, old_ode.u_vars) ≡ map(string, new_ode.u_vars)
        # @test map(string, old_ode.y_vars) ≡ map(string, new_ode.y_vars)
        @test length(old_ode.y_vars) == length(new_ode.y_vars)
        @test length(old_ode.u_vars) == length(new_ode.u_vars)
        @test base_ring(parent(first(values(var_mapping)))) == parent(old_ode)

        # We check that the PS solutions of the new_ode coincide with the solutions
        # of the old_ode projected onto the new variables
        ord = 5
        # NOTE: there may be an unlucky specialization point, and newton iteration
        # will fail, and tests will fail as a result
        bound = 100
        old_params = old_ode.parameters
        old_vars = old_ode.x_vars
        old_inputs = old_ode.u_vars
        param_spec_point = map(Nemo.QQ, rand(1:bound, length(old_params)))
        old_param_spec = Dict(old_params .=> param_spec_point)
        ic_point = map(Nemo.QQ, rand(1:bound, length(old_vars)))
        old_var_ic = Dict(old_vars .=> ic_point)
        input_ts = map(_ -> [Nemo.QQ(rand(1:bound))], 1:length(old_inputs))
        old_input_ts =
            Dict{eltype(old_vars), Vector{Nemo.QQFieldElem}}(old_inputs .=> input_ts)

        old_solutions = StructuralIdentifiability.power_series_solution(
            old_ode,
            old_param_spec,
            old_var_ic,
            old_input_ts,
            ord,
        )

        new_params = new_ode.parameters
        new_vars = new_ode.x_vars
        new_inputs = new_ode.u_vars
        new_param_spec = Dict(
            new_param => StructuralIdentifiability.eval_at_dict(
                var_mapping[new_param],
                old_param_spec,
            ) for new_param in new_params
        )
        new_var_ic = Dict(
            new_var => StructuralIdentifiability.eval_at_dict(
                var_mapping[new_var],
                merge(old_param_spec, old_var_ic),
            ) for new_var in new_vars
        )
        new_input_ts = Dict{eltype(new_vars), Vector{Nemo.QQFieldElem}}(
            new_input => old_input_ts[numerator(var_mapping[new_input])] for
            new_input in new_inputs
        )

        new_solutions = StructuralIdentifiability.power_series_solution(
            new_ode,
            new_param_spec,
            new_var_ic,
            new_input_ts,
            ord,
        )

        # NOTE: test that y's are mapped one to one!
        @test map(string, [var_mapping[output] for output in new_ode.y_vars]) ≡ map(string, old_ode.y_vars)

        # Test IO behavior 
        for var in vcat(new_ode.y_vars)
            new_solution = new_solutions[var]
            old_solution = old_solutions[numerator(var_mapping[var])]
            @test new_solution == old_solution
        end

        # Test state dynamics
        projected_solutions = Dict(
            var =>
                StructuralIdentifiability.eval_at_dict(var_mapping[var], old_solutions)
            for var in new_ode.x_vars
        )
        for var in vcat(new_ode.x_vars)
            new_solution = new_solutions[var]
            prj_solution = projected_solutions[var]
            @test new_solution == prj_solution
        end

        # Test relations
        for relation in implicit_relations
            @test iszero(StructuralIdentifiability.eval_at_dict(relation, var_mapping))
        end

        return nothing
    end

    cases = [
        (
            ode = StructuralIdentifiability.@ODEmodel(
                x1'(t) = x1(t) + a + u1(t),
                x2'(t) = x2(t) + b + u2(t),
                y1(t) = x1(t),
                y2(t) = x2(t),
            ),
        ),
        (
            ode = StructuralIdentifiability.@ODEmodel(
                x1'(t) = a * x1,
                x2'(t) = b * x2,
                y(t) = x1 + x2
            ),
        ),
        (
            ode = StructuralIdentifiability.@ODEmodel(
                x1'(t) = (a + b) * x1 // x2 + a,
                x2'(t) = x1 // (a + c) + c,
                y(t) = x1
            ),
        ),
        (
            ode = StructuralIdentifiability.@ODEmodel(
                T1'(t) = a * T1 - b * T1 * T2 + u(t),
                T2'(t) = -c * T2 + d * T1 * T2,
                y(t) = T1
            ),
        ),
        (
            #=
            Identifiable functions:
                x1, Θ * x2

            Under
                T1 = x1, 
                T2 = Θ * x2

            Reparametrizes into:
                T1' = T1 + T2,
                T2' = 0,
                y   = T1
            =#
            ode = StructuralIdentifiability.@ODEmodel(
                x1'(t) = x1(t) + Θ * x2(t),
                x2'(t) = 0,
                y(t) = x1(t)
            ),
        ),
        (
            #=
            Identifiable functions:
                x1, α * x2, α * x3, α * C

            Under
                T1 = α * x2, 
                T2 = α * x3,
                T3 = x1,
                T4 = α * C

            Reparametrizes into:
                T3'(t) = T1(t)
                T2'(t) = T4
                T1'(t) = T2(t)
                y(t)   = T3(t)
            =#
            ode = StructuralIdentifiability.@ODEmodel(
                x1'(t) = α * x2(t),
                x2'(t) = x3(t),
                x3'(t) = C,
                y(t) = x1(t)
            ),
        ),
        (
            #=
            Identifiable functions:
                x1^2 + x2^2, α

            Under
                T1 = x1^2 + x2^2,
                T2 = α

            Reparametrizes into:
                T1' = 2*T1*T2
                y   = 1//2*T1
            =#
            ode = StructuralIdentifiability.@ODEmodel(
                x1'(t) = α * (x1 - x2),
                x2'(t) = α * (x1 + x2),
                y(t) = (x1^2 + x2^2) // 2,
            ),
        ),
        (
            # Pivastatin
            # Reparametrizes into a 4-dimensional nonlinear model with no algebraic
            # relations -- how is this even legal?? 
            ode = StructuralIdentifiability.@ODEmodel(
                x1'(t) =
                    k3 * x3(t) - r3 * x1(t) - k1 * x1(t) * (T0 - x2(t)) + r1 * x2(t),
                x2'(t) = k1 * x1(t) * (T0 - x2(t)) - (r1 + k2) * x2(t),
                x3'(t) = r3 * x1(t) - (k3 + k4) * x3(t) + k2 * x2(t),
                y1(t) = k * (x2(t) + x3(t))
            ),
        ),
        (
            # Transfection_4State
            ode = StructuralIdentifiability.@ODEmodel(
                mRNA'(t) = -d1 * mRNA(t) - d2 * mRNA(t) * enz(t),
                GFP'(t) = kTL * mRNA(t) - b * GFP(t),
                enz'(t) = d3 * mRNAenz(t) - d2 * mRNA(t) * enz(t),
                mRNAenz'(t) = -d3 * mRNAenz(t) + d2 * mRNA(t) * enz(t),
                y1(t) = GFP(t)
            ),
        ),
        (
            # LLW1987_io
            ode = StructuralIdentifiability.@ODEmodel(
                x1'(t) = -p1 * x1(t) + p2 * u(t),
                x2'(t) = -p3 * x2(t) + p4 * u(t),
                x3'(t) = -(p1 + p3) * x3(t) + (p4 * x1(t) + p2 * x2(t)) * u(t),
                y1(t) = x3(t)
            ),
        ),
        (
            # Take care of monomial orders in GB!
            ode = StructuralIdentifiability.@ODEmodel(
                x1'(t) = x1 + x2^2 + a^2,
                x2'(t) = x2 + a * d^3,
                y(t) = x1
            ),
        ),
        (
            # SLIQR system
            ode = StructuralIdentifiability.@ODEmodel(
                S'(t) = -b * In(t) * S(t) * Ninv - u(t) * S(t) * Ninv,
                L'(t) = b * In(t) * S(t) * Ninv - a * L(t),
                In'(t) = a * L(t) - g * In(t) + s * Q(t),
                Q'(t) = (1 - e) * g * In(t) - s * Q(t),
                y(t) = In(t) * Ninv
            ),
        ),
        (
            ode = StructuralIdentifiability.@ODEmodel(
                x1'(t) = a * x1 + b * x2 + u(t),
                x2'(t) = b * x1 + c * x2,
                y(t) = x1
            ),
        ),
    ]

    @testset "Global reparametrizations" begin
        # Test that variables are mapped properly
        ode = cases[1].ode
        (new_ode, new_vars, implicit_relations) =
            StructuralIdentifiability.reparametrize_global(ode)
        @test length(new_vars) == length(gens(parent(new_ode))) == 8
        @test length(new_ode.x_vars) == 2
        # @test map(string, new_ode.y_vars) ≡ map(string, ode.y_vars)
        # @test map(string, new_ode.u_vars) ≡ map(string, ode.u_vars)
        @test length(new_ode.parameters) == 2
        @test isempty(implicit_relations)

        # Test that relations are functional
        ode = cases[2].ode
        (new_ode, new_vars, implicit_relations) =
            StructuralIdentifiability.reparametrize_global(ode)
        @test length(implicit_relations) == 1

        for case in cases
            ode = case.ode
            (new_ode, new_vars, implicit_relations) =
                StructuralIdentifiability.reparametrize_global(ode)
            test_reparametrization(ode, new_ode, new_vars, implicit_relations)
        end
    end
end
