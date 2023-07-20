@testset "Generators of observable states" begin
    ode = @ODEmodel(
        x1'(t) = a * x1(t) - b * x1(t) * x2(t),
        x2'(t) = -c * x2(t) + d * x1(t) * x2(t),
        y(t) = x1(t)
    )

    sg = states_generators(ode, find_ioequations(ode))
    @test sg == [
        x1 // 1,
        (a * x1 - b * x1 * x2) // 1,
        ((a - b * x2) * (a * x1 - b * x1 * x2) - b * x1 * (-c * x2 + d * x1 * x2)) // 1,
    ]

    ode = @ODEmodel(
        x1'(t) = x2(t),
        x2'(t) = x3(t),
        x3'(t) = x4(t),
        x4'(t) = 1 // x2(t),
        y(t) = x1(t)
    )
    sg = states_generators(ode, find_ioequations(ode))
    @test sg == [x1 // 1, x2 // 1, x3 // 1, x4 // 1, x2 // 1]

    ode = @ODEmodel(
        x1'(t) = x2(t),
        x2'(t) = (x1(t) + a * x2(t) * u(t)) // (u(t)^2 + x1(t)),
        y(t) = x1(t)
    )
    sg = states_generators(ode, find_ioequations(ode))
    @test Set(sg) == Set([x1 // 1, x2 // 1, x1 // 1, a * x2 // 1, x1 // 1])
end
