# using Profile
using IterTools
using Oscar

function rand_poly(deg, vars)
    result = 0
    indices = vcat(collect(1:length(vars)), collect(1:length(vars)))
    monomials = []
    for d in 0:deg
        for subs in IterTools.subsets(indices, d)
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
    exp_t = ps_ode_solution([x_dot - x], Dict{fmpq_mpoly, fmpq}(x => 1), Dict{fmpq_mpoly, Array{fmpq, 1}}(), 20)[x]
    @test valuation(ps_diff(exp_t) - exp_t) == 19

    R, (x, y, x_dot, y_dot, u) = Nemo.PolynomialRing(Nemo.QQ, ["x", "y", "x_dot", "y_dot", "u"])
    prec = 100
    eqs = [
        x_dot - x + 3 * x * y - u,
        y_dot + 2 * y - 4 * x * y
    ]
    u_coeff = [rand(1:5) for i in 1:prec]
    sol = ps_ode_solution(eqs, Dict(x => 0, y => -2), Dict(u => u_coeff), prec)
    @test map(e -> valuation(evaluate(e, [sol[v] for v in gens(R)])), eqs) == [prec - 1, prec - 1]

    # Testing the function ps_ode_solution by itself
    for i in 1:3
        prec = 100
        varnames = vcat(
            ["x_$(i)_dot" for i in 1:5],
            ["x_$i" for i in 1:5],
            ["u_$i" for i in 1:3],
        )
        R, vars = Nemo.PolynomialRing(Nemo.GF(2^31 - 1), varnames)
        eqs = [rand_poly(1, vars[6:end]) * vars[i] - rand_poly(2, vars[6:end]) for i in 1:5]
        ic = Dict(vars[i + 5] => rand(-5:5) for i in 1:5)
        inputs = Dict(u => [rand(-3:3) for i in 1:prec] for u in vars[11:end])
        @time sol = ps_ode_solution(eqs, ic, inputs, prec)
        evals = map(e -> valuation(evaluate(e, [sol[v] for v in gens(R)])), eqs)
        for e in evals
            @test e >= prec - 1
        end
    end

    # Testing ps_ode_solution in conjuntion with the ODE class
    for i in 1:3
        prec = 100
        varnames = vcat(
            ["x_$i" for i in 1:3],
            ["p_$i" for i in 1:3],
            ["u_$i" for i in 1:2],
        )
        R, vars = Nemo.PolynomialRing(Nemo.GF(2^31 - 1), varnames)
        PType = gfp_mpoly
        TDict = Dict{PType, Union{PType, Generic.Frac{PType}}}
        eqs = TDict(vars[i] => rand_poly(2, vars) // rand_poly(1, vars) for i in 1:3)
        ode = ODE{PType}(eqs, TDict(), vars[(end - 1):end])
        ic = Dict(vars[i] => rand(-5:5) for i in 1:3)
        param_vals = Dict(vars[i + 3] => rand(-5:5) for i in 1:3)
        inputs = Dict(u => [rand(-3:3) for i in 1:prec] for u in vars[7:end])
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

#Profile.print(format=:flat, maxdepth=4, sortedby=:count)
end
