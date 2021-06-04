@testset "Generating polynomial system via Lie derivatives" begin
 
    ode = @ODEmodel(
        x'(t) = a * x(t) + b,
        y(t) = c * x(t)
    )

    ps_sol = power_series_solution(ode, Dict(a => 3, b => 5, c => -2), Dict(x => 10), Dict{fmpq_mpoly, Array{Int, 1}}(), 6)
    ps_sol = Dict(var_to_str(v) => sol for (v, sol) in ps_sol)

    result = get_poly_system_Lie(ode)

    @test length(gens(result[:diff_ring].ring)) == 9

    a, b, c, x, y0, y1, y2, y3, y4 = map(
        v -> str_to_var(v, result[:diff_ring].ring),
        vcat(["a", "b", "c", diffvar("x", 0)], [diffvar("y", i) for i in 0:4])
    )

    @test result[:square_system] == [
        y0 - c * x,
        y1 - a * c * x - b * c,
        y2 - a^2 * c * x - a * b * c,
        y3 - a^3 * c * x - a^2 * b * c
    ]

    @test result[:extra_equations] == [y4 - a^4 * c * x - a^3 * b * c]

    eval_point = produce_point(ps_sol, result[:diff_ring])
    eval_point = [eval_point[x] for x in gens(result[:diff_ring].ring)]
    for eq in vcat(result[:square_system], result[:extra_equations])
        @test evaluate(eq, eval_point) == 0
    end
    

    #----------------------------

    result_lazy = get_poly_system_lazy(ode)
    @test length(gens(result_lazy[:diff_ring].ring)) == 13

    eval_point = produce_point(ps_sol, result_lazy[:diff_ring])
    eval_point = [eval_point[x] for x in gens(result_lazy[:diff_ring].ring)]
    for eq in vcat(result_lazy[:square_system], result_lazy[:extra_equations])
        @test evaluate(eq, eval_point) == 0
    end
 
end
