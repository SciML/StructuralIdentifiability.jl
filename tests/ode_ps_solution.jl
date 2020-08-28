# using Profile

function rand_poly(deg, vars)
    result = 0
    iters = [0:deg for v in vars]
    for exp in Iterators.product(iters...)
        monomial = rand(-5:5)
        for i in 1:length(vars)
            monomial *= vars[i]^exp[i]
        end
        result += monomial
    end
    return result
end

@testset "Power series solution for an ODE system" begin

    R, (x, x_dot) = PolynomialRing(QQ, ["x", "x_dot"])
    exp_t = ps_ode_solution([x_dot - x], Dict(x => 1), Dict(), 20)[x]
    @test valuation(ps_diff(exp_t) - exp_t) == 19

    R, (x, y, x_dot, y_dot, u) = PolynomialRing(QQ, ["x", "y", "x_dot", "y_dot", "u"])
    prec = 100
    eqs = [
        x_dot - x + 3 * x * y - u,
        y_dot + 2 * y - 4 * x * y
    ]
    u_coeff = [rand(1:5) for i in 1:prec]
    sol = ps_ode_solution(eqs, Dict(x => 0, y => -2), Dict(u => u_coeff), prec)
    @test map(e -> valuation(evaluate(e, [sol[v] for v in gens(R)])), eqs) == [prec - 1, prec - 1]

    for i in 1:3
        prec = 100
        varnames = vcat(
            ["x_$(i)_dot" for i in 1:5],
            ["x_$i" for i in 1:5],
            ["u_$i" for i in 1:3],
        )
        R, vars = PolynomialRing(GF(2^31 - 1), varnames)
        eqs = [rand_poly(1, vars[6:end]) * vars[i] - rand_poly(2, vars[6:end]) for i in 1:5]
        ic = Dict(vars[i + 5] => rand(-5:5) for i in 1:5)
        inputs = Dict(u => [rand(-3:3) for i in 1:prec] for u in vars[11:end])
        @time sol = ps_ode_solution(eqs, ic, inputs, prec)
        evals = map(e -> valuation(evaluate(e, [sol[v] for v in gens(R)])), eqs)
        for e in evals
            @test e >= prec - 1
        end
    end

#Profile.print(format=:flat, maxdepth=4, sortedby=:count)
end
