include(joinpath(@__DIR__, "..", "shared", "test_setup.jl"))

@testset "Computing extra projection" begin

    ode = @ODEmodel(
        x1'(t) = (1 + x1(t)^2) // 2,
        x2'(t) = (1 - x1(t)^2) // (1 + x1(t)^2),
        y1(t) = 2 * x1(t) // (b * (1 + x1(t)^2)),
        y2(t) = x2(t)
    )
    ioeq_dict = find_ioequations(ode)
    ioeq_dict = find_extra_projection(ode, ioeq_dict)
    @test length(ioeq_dict) == 3
    @test [degree(eq, v) for (v, eq) in ioeq_dict] == [2, 2, 2]
    @test Set([degree(eq, str_to_var("b", parent(eq))) for eq in values(ioeq_dict)]) == Set(2, 4)
end
