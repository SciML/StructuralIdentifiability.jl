using Test
using TestSetExtensions

@testset "Lie derivative" begin

    ode = @ODEmodel(
        x1'(t) = a * x1(t) + b * u(t),
        x2'(t) = 1 // (1 - x1(t) - x2(t)),
        y(t) = x1(t)
    )

    @test lie_derivative(x1 + x2, ode) == ((a * x1 + b * u) * (1 - x1 - x2) + 1) // (1 - x1 - x2)

    ode = @ODEmodel(
        x1'(t) = x2(t) // x1(t) + x2(t) * u(t),
        x2'(t) = -1 - x1(t) * u(t),
        y(t) = 1
    )

    @test lie_derivative(x1^2 + x2^2 // 1, ode) == zero(parent(ode)) // 1

end
