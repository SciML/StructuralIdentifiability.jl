if GROUP == "All" || GROUP == "Core"
    #------- Auxiliary functions --------------------------------------------------

    function diff_sol_Lie_derivatives(ode::ODE, params, ic, inputs, prec::Int)
        # creating a new ring with variables for the derivatives of u
        new_varnames = map(var_to_str, gens(ode.poly_ring))
        if length(ode.u_vars) > 0
            append!(
                new_varnames,
                [var_to_str(u) * "_$i" for u in ode.u_vars for i in 0:(prec - 1)],
            )
        end
        new_ring, vars = Nemo.polynomial_ring(base_ring(ode.poly_ring), new_varnames)

        # mapping everything to the new ring
        eval_point = Dict(v => switch_ring(v, new_ring) for v in gens(ode.poly_ring))
        for u in ode.u_vars
            eval_point[u] = str_to_var(var_to_str(u) * "_0", new_ring)
        end

        new_eqs = Dict()
        for (x, f) in ode.x_equations
            new_eqs[str_to_var(var_to_str(x), new_ring)] = eval_at_dict(f, eval_point)
        end
        params, ic = map(
            d -> Dict(str_to_var(string(k), new_ring) => v for (k, v) in d),
            [params, ic],
        )

        # computing Lie derivatives
        derivation = copy(new_eqs)
        for u in ode.u_vars
            for i in 0:(prec - 2)
                derivation[str_to_var(var_to_str(u) * "_$i", new_ring)] =
                    str_to_var(var_to_str(u) * "_$(i + 1)", new_ring)
            end
        end
        Lie_derivatives = Dict()
        for (y, g) in ode.y_equations
            Lie_derivatives[y] = Array{Any, 1}([eval_at_dict(g, eval_point)])
            for i in 1:prec
                push!(
                    Lie_derivatives[y],
                    sum(
                        derivative(Lie_derivatives[y][end], v) * get(derivation, v, 0) for
                        v in gens(new_ring)
                    ),
                )
            end
        end

        # producing the result
        eval_point = merge(params, ic)
        for u in ode.u_vars
            for i in 1:prec
                eval_point[str_to_var(var_to_str(u) * "_$(i - 1)", new_ring)] =
                    inputs[u][i] * factorial(i - 1)
            end
        end

        result = Dict()
        for y in ode.y_vars
            result[y] = Dict()
            for v in vcat(ode.x_vars, ode.parameters)
                result[y][v] = []
                for j in 1:prec
                    push!(
                        result[y][v],
                        eval_at_dict(
                            derivative(Lie_derivatives[y][j], str_to_var("$v", new_ring)),
                            eval_point,
                        ),
                    )
                end
            end
        end

        return result
    end

    #------------------------------------------------------------------------------

    function rand_poly(deg, vars)
        if deg == 0
            return parent(vars[1])(1)
        end
        result = 0
        indices = collect(1:length(vars))
        monomials = []
        for d in 0:deg
            for subs in StructuralIdentifiability.IterTools.subsets(indices, d)
                push!(monomials, subs)
            end
        end

        for subs in monomials
            monom = rand(-50:50)
            for v_ind in subs
                monom *= vars[v_ind]
            end
            result += monom
        end

        return result
    end

    #------------------------------------------------------------------------------

    @testset "Partial derivatives of an output w.r.t. to initial conditions and parameters" begin
        test_cases = []
        P = QQMPolyRingElem
        DType = Union{P, Generic.Frac{P}}

        ode = @ODEmodel(x'(t) = x(t) + a, y(t) = x(t)^2)
        push!(
            test_cases,
            Dict(
                :ODE => ode,
                :ic => Dict(x => Nemo.QQ(rand(1:10))),
                :param_vals => Dict(a => Nemo.QQ(rand(1:10))),
                :inputs => Dict{P, Array{QQFieldElem, 1}}(),
                :prec => 20,
            ),
        )

        ode = @ODEmodel(x'(t) = x(t)^2 + a, y1(t) = x(t) + a^2, y2(t) = x(t)^3)
        push!(
            test_cases,
            Dict(
                :ODE => ode,
                :ic => Dict(x => Nemo.QQ(rand(1:10))),
                :param_vals => Dict(a => Nemo.QQ(rand(1:10))),
                :inputs => Dict{P, Array{QQFieldElem, 1}}(),
                :prec => 20,
            ),
        )

        ode = @ODEmodel(
            x'(t) = x(t)^2 + 2 * x(t) * y(t) - 3 * a * y(t),
            y'(t) = x(t)^2 + a * b - b^2 + 4 * b * x(t),
            y1(t) = a * x(t),
            y2(t) = b * y(t)^2 - y(t)
        )
        push!(
            test_cases,
            Dict(
                :ODE => ode,
                :ic => Dict(x => Nemo.QQ(rand(1:10)), y => Nemo.QQ(rand(1:10))),
                :param_vals => Dict(a => Nemo.QQ(rand(1:10)), b => Nemo.QQ(rand(1:10))),
                :inputs => Dict{P, Array{QQFieldElem, 1}}(),
                :prec => 8,
            ),
        )

        ode = @ODEmodel(x'(t) = u(t) + a, y(t) = x(t))
        push!(
            test_cases,
            Dict(
                :ODE => ode,
                :ic => Dict(x => Nemo.QQ(rand(1:10))),
                :param_vals => Dict(a => Nemo.QQ(rand(1:10))),
                :inputs => Dict(u => [Nemo.QQ(rand(-3:3)) for i in 1:20]),
                :prec => 20,
            ),
        )

        F = Nemo.Native.GF(2^31 - 1)
        P = gfp_mpoly
        DType = Union{P, Generic.Frac{P}}

        varnames = vcat(
            ["x_$i" for i in 1:3],
            ["p_$i" for i in 1:3],
            ["u_$i" for i in 1:2],
            ["y_$i" for i in 1:3],
        )
        R, vars = Nemo.polynomial_ring(F, varnames)
        push!(
            test_cases,
            Dict(
                :ODE => ODE{P}(
                    Dict{P, DType}(vars[i] => rand_poly(1, vars[1:8]) for i in 1:3),
                    Dict{P, DType}(vars[i] => rand_poly(2, vars[1:8]) for i in 9:11),
                    vars[7:8],
                ),
                :ic => Dict(vars[i] => F(rand(1:50)) for i in 1:3),
                :param_vals => Dict(vars[i + 3] => F(rand(1:50)) for i in 1:3),
                :inputs => Dict(u => [F(rand(-30:30)) for i in 1:6] for u in vars[7:8]),
                :prec => 6,
            ),
        )

        varnames = vcat(
            ["x_$i" for i in 1:3],
            ["p_$i" for i in 1:3],
            ["u_$i" for i in 1:2],
            ["y_$i" for i in 1:3],
        )
        R, vars = Nemo.polynomial_ring(F, varnames)
        push!(
            test_cases,
            Dict(
                :ODE => ODE{P}(
                    Dict{P, DType}(vars[i] => rand_poly(2, vars[1:8]) for i in 1:3),
                    Dict{P, DType}(vars[i] => rand_poly(2, vars[1:8]) for i in 9:11),
                    vars[7:8],
                ),
                :ic => Dict(vars[i] => F(rand(1:50)) for i in 1:3),
                :param_vals => Dict(vars[i + 3] => F(rand(1:50)) for i in 1:3),
                :inputs => Dict(u => [F(rand(-30:30)) for i in 1:6] for u in vars[7:8]),
                :prec => 6,
            ),
        )

        varnames = vcat(["x_$i" for i in 1:2], ["p_$i" for i in 1:2], "u", ["y_1", "y_2"])
        R, vars = Nemo.polynomial_ring(F, varnames)
        push!(
            test_cases,
            Dict(
                :ODE => ODE{P}(
                    Dict{P, DType}(
                        vars[i] => rand_poly(1, vars[1:5]) // (vars[1] + vars[3]) for
                        i in 1:2
                    ),
                    Dict{P, DType}(vars[i] => rand_poly(1, vars[1:5]) for i in 6:7),
                    [vars[5]],
                ),
                :ic => Dict(vars[i] => F(rand(1:50)) for i in 1:2),
                :param_vals => Dict(vars[i + 2] => F(rand(1:50)) for i in 1:2),
                :inputs => Dict(vars[5] => [F(rand(-30:30)) for i in 1:4]),
                :prec => 4,
            ),
        )

        for case in test_cases
            ode, prec = case[:ODE], case[:prec]
            @time sol1 =
                differentiate_output(ode, case[:param_vals], case[:ic], case[:inputs], prec)
            sol2 = diff_sol_Lie_derivatives(
                ode,
                case[:param_vals],
                case[:ic],
                case[:inputs],
                prec,
            )
            for y in ode.y_vars
                for v in vcat(ode.x_vars, ode.parameters)
                    @test sol2[y][v] == [
                        base_ring(ode.poly_ring)(coeff(sol1[y][v], j) * factorial(j)) for
                        j in 0:(prec - 1)
                    ]
                end
            end
        end
    end
end
