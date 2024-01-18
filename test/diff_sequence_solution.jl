if GROUP == "All" || GROUP == "Core"
    @testset "Computing variations around a sequence solution" begin
        # Computing sensitivities directly be explicitly writing down Lie derivatives
        function use_lie_derivatives(
            dds::ODE{P},
            params::Dict{P, T},
            ic::Dict{P, T},
            inputs::Dict{P, Array{T, 1}},
            num_terms::Int,
        ) where {T <: Generic.FieldElem, P <: MPolyElem{T}}
            newvars = [var_to_str(v) for v in gens(dds.poly_ring)]
            append!(
                newvars,
                [var_to_str(v) * "$i" for v in dds.u_vars for i in 1:num_terms],
            )
            R, _ = StructuralIdentifiability.Nemo.PolynomialRing(
                base_ring(dds.poly_ring),
                newvars,
            )
            explicit_sol = merge(
                Dict(
                    parent_ring_change(x, R) => Vector{Any}([parent_ring_change(x, R)])
                    for (x, eq) in dds.x_equations
                ),
                Dict(
                    parent_ring_change(y, R) =>
                        Vector{Any}([parent_ring_change(eq, R)]) for
                    (y, eq) in dds.y_equations
                ),
            )
            time_step = merge(
                Dict(
                    parent_ring_change(x, R) => parent_ring_change(eq, R) for
                    (x, eq) in dds.x_equations
                ),
                Dict(
                    parent_ring_change(p, R) => parent_ring_change(p, R) for
                    p in dds.parameters
                ),
                Dict(
                    parent_ring_change(u, R) => str_to_var(var_to_str(u) * "1", R) for
                    u in dds.u_vars
                ),
                Dict(
                    str_to_var(s * "$i", R) => str_to_var(s * "$(i + 1)", R) for
                    s in map(var_to_str, dds.u_vars) for i in 1:(num_terms - 1)
                ),
            )
            eval_dict = merge(
                Dict(parent_ring_change(p, R) => v for (p, v) in params),
                Dict(parent_ring_change(x, R) => v for (x, v) in ic),
                Dict(parent_ring_change(u, R) => val[1] for (u, val) in inputs),
                Dict(
                    str_to_var(var_to_str(u) * "$i", R) => inputs[u][i + 1] for
                    u in dds.u_vars for i in 1:(num_terms - 1)
                ),
            )
            generalized_parameters =
                [parent_ring_change(p, R) for p in vcat(dds.x_vars, dds.parameters)]
            for i in 2:num_terms
                for (k, v) in explicit_sol
                    push!(explicit_sol[k], eval_at_dict(v[end], time_step))
                end
            end
            part_diffs = Dict(
                (f, p) => [] for f in keys(explicit_sol) for p in generalized_parameters
            )
            for i in 1:num_terms
                for (f, ders) in explicit_sol
                    for p in generalized_parameters
                        push!(
                            part_diffs[(f, p)],
                            eval_at_dict(derivative(ders[i], p), eval_dict),
                        )
                    end
                end
            end
            return Dict(
                (
                    parent_ring_change(k[1], dds.poly_ring),
                    parent_ring_change(k[2], dds.poly_ring),
                ) => res for (k, res) in part_diffs
            )
        end

        locQQ = StructuralIdentifiability.Nemo.QQ

        ode = @ODEmodel(a'(t) = a(t)^2 + b, y(t) = 1 / (a(t) * c(t)))
        params = Dict(b => locQQ(1))
        ic = Dict(a => locQQ(2))
        inputs = Dict(c => [locQQ(1), locQQ(-2), locQQ(3), locQQ(-4), locQQ(5)])
        seq_sol, diff_sol = differentiate_sequence_solution(ode, params, ic, inputs, 4)
        diff_y = differentiate_sequence_output(ode, params, ic, inputs, 4)
        lie_ders_sol = use_lie_derivatives(ode, params, ic, inputs, 4)
        merged = merge(diff_sol, diff_y)
        @test merged == lie_ders_sol

        ode = @ODEmodel(
            a'(t) =
                (23 * k1 * a(t) - 7 * b(t)^3) // (a(t)^2 + b(t)^2) - c(t)^3 * k1 * b(t),
            b'(t) = a(t) + 17 * (b(t) - c(t))^2 + 1 // (a(t) + b(t) - k2),
            y(t) = (a(t) + b(t) - c(t)) // (a(t)^2 + k2)
        )
        params = Dict(k1 => locQQ(1), k2 => locQQ(2))
        ic = Dict(a => locQQ(3), b => locQQ(-4))
        inputs = Dict(c => [locQQ(5), locQQ(-6), locQQ(7), locQQ(-8)])
        seq_sol, diff_sol = differentiate_sequence_solution(ode, params, ic, inputs, 2)
        diff_y = differentiate_sequence_output(ode, params, ic, inputs, 2)
        lie_ders_sol = use_lie_derivatives(ode, params, ic, inputs, 2)
        merged = merge(diff_sol, diff_y)
        @test merged == lie_ders_sol
    end
end
