include(joinpath(@__DIR__, "..", "shared", "test_setup.jl"))

@testset "Identifiable functions of parameters" begin
    ode = @ODEmodel(
        x1'(t) = -(k21 + V_M / (K_M + x1(t))) * x1(t) + k12 * x2(t) + b1 * u(t),
        x2'(t) = k21 * x1(t) - (k02 + k12) * x2(t),
        y(t) = c * x1(t)
    )

    @test observation_field(ode) == find_identifiable_functions(ode, with_states = true)
end
