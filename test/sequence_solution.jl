if GROUP == "All" || GROUP == "Core"
    @testset "Sequence solutions in the discrete case" begin
        # Fibonacci example
        ode = @DDSmodel(f1(t + 1) = f1(t) + f0(t), f0(t + 1) = f1(t), y(t) = f1(t))

        sol = sequence_solution(
            ode,
            Dict{QQMPolyRingElem, Int}(),
            Dict(f0 => 1, f1 => 1),
            Dict{QQMPolyRingElem, Array{Int, 1}}(),
            10,
        )
        @test sol[f1] == [1, 2, 3, 5, 8, 13, 21, 34, 55, 89]
        @test sol[f0] == [1, 1, 2, 3, 5, 8, 13, 21, 34, 55]

        # Sum of powers
        ode = @DDSmodel(a(t + 1) = 2 * a(t) + b * u(t), y(t) = a(t))
        sol = sequence_solution(
            ode,
            Dict(b => 3),
            Dict(a => 1),
            Dict(u => [3^i for i in 0:9]),
            10,
        )
        @test sol[a] == [-2^i + 3^i for i in 1:10]
    end
end
