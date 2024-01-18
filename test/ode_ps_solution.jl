if GROUP == "All" || GROUP == "Core"
    function rand_poly(deg, vars)
        result = 0
        indices = vcat(collect(1:length(vars)), collect(1:length(vars)))
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

    @testset "Power series solution for an ODE system" begin
        R, (x, x_dot) = Nemo.PolynomialRing(Nemo.QQ, ["x", "x_dot"])
        exp_t = ps_ode_solution(
            [x_dot - x],
            Dict{fmpq_mpoly, fmpq}(x => 1),
            Dict{fmpq_mpoly, Array{fmpq, 1}}(),
            20,
        )[x]
        @test valuation(ps_diff(exp_t) - exp_t) == 19

        R, (x, y, x_dot, y_dot, u) =
            Nemo.PolynomialRing(Nemo.QQ, ["x", "y", "x_dot", "y_dot", "u"])
        prec = 100
        eqs = [x_dot - x + 3 * x * y - u, y_dot + 2 * y - 4 * x * y]
        u_coeff = [rand(1:5) for i in 1:prec]
        sol = ps_ode_solution(eqs, Dict(x => 0, y => -2), Dict(u => u_coeff), prec)
        @test map(e -> valuation(evaluate(e, [sol[v] for v in gens(R)])), eqs) == [prec - 1, prec - 1]

        F = Nemo.Native.GF(2^31 - 1)
        prec = 100

        # Testing the function ps_ode_solution by itself
        for i in 1:30
            # Setting up the ring
            NUMX = 5
            NUMU = 3
            varnames = vcat(
                ["x_$(i)_dot" for i in 1:NUMX],
                ["x_$i" for i in 1:NUMX],
                ["u_$i" for i in 1:NUMU],
            )
            R, vars = Nemo.PolynomialRing(F, varnames)

            # Generating the initial conditions and inputs
            ic = Dict(vars[i + NUMX] => F(rand(-5:5)) for i in 1:NUMX)
            inputs =
                Dict(u => [F(rand(-3:3)) for i in 1:prec] for u in vars[(2 * NUMX + 1):end])
            # a dictionary for evaluation at zero (to avoid singularities)
            eval_at_zero = merge(ic, Dict(u => val[1] for (u, val) in inputs))

            # Generating denominators not vanishing at t = 0
            denominators = [rand_poly(1, vars[(NUMX + 1):end]) for i in 1:NUMX]
            while any([eval_at_dict(p, eval_at_zero) == 0 for p in denominators])
                denominators = [rand_poly(1, vars[(NUMX + 1):end]) for i in 1:NUMX]
            end

            eqs = [
                denominators[i] * vars[i] - rand_poly(2, vars[(NUMX + 1):end]) for
                i in 1:NUMX
            ]
            @time sol = ps_ode_solution(eqs, ic, inputs, prec)
            evals = map(e -> valuation(evaluate(e, [sol[v] for v in gens(R)])), eqs)
            for e in evals
                @test e >= prec - 1
            end
        end

        # Testing ps_ode_solution in conjuntion with the ODE class
        for i in 1:30
            # Setting up the ring
            NUMX = 3
            NUMP = 3
            NUMU = 2
            varnames = vcat(
                ["x_$i" for i in 1:NUMX],
                ["p_$i" for i in 1:NUMP],
                ["u_$i" for i in 1:NUMU],
            )
            R, vars = Nemo.PolynomialRing(F, varnames)
            PType = gfp_mpoly
            TDict = Dict{PType, Union{PType, Generic.Frac{PType}}}

            # Generating the intial conditions etc
            ic = Dict(vars[i] => F(rand(-5:5)) for i in 1:NUMX)
            param_vals = Dict(vars[i + NUMX] => F(rand(-5:5)) for i in 1:NUMP)
            inputs = Dict(
                u => [F(rand(-3:3)) for i in 1:prec] for u in vars[(NUMX + NUMP + 1):end]
            )
            eval_at_zero = merge(ic, param_vals, Dict(u => val[1] for (u, val) in inputs))

            # Generating denominators not vanishing at t = 0
            denominators = [rand_poly(1, vars) for i in 1:NUMX]
            while any([eval_at_dict(p, eval_at_zero) == 0 for p in denominators])
                denominators = [rand_poly(1, vars) for i in 1:NUMX]
            end

            eqs = TDict(vars[i] => rand_poly(2, vars) // denominators[i] for i in 1:NUMX)
            ode = ODE{PType}(eqs, TDict(), vars[(NUMX + NUMP + 1):end])
            @time sol = power_series_solution(ode, param_vals, ic, inputs, prec)
            ps_ring = parent(collect(values(sol))[1])
            for (p, val) in param_vals
                sol[p] = ps_ring(val)
            end
            eval_point = [sol[v] for v in vars]
            for (xd, f) in eqs
                num, den = unpack_fraction(f)
                lhs = divexact(evaluate(num, eval_point), evaluate(den, eval_point))
                rhs = ps_diff(sol[xd])
                @test valuation(lhs - rhs) >= prec - 1
            end
        end
    end
end
