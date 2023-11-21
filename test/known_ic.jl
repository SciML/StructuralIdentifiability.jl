cases = []

ode = @ODEmodel(
    x1'(t) = a * x1(t) - b * x1(t) * x2(t),
    x2'(t) = -c * x2(t) + d * x1(t) * x2(t),
    y(t) = x1(t)
)

push!(cases, Dict(
    :ode => ode,
    :known => [x2],
    :correct => [a, b, c, d, x1, x2],
))

ode = @ODEmodel(
    x1'(t) = a  + x2(t) + x3(t),
    x2'(t) = b,
    x3'(t) = c,
    y(t) = x1(t)
)

push!(cases, Dict(
    :ode => ode,
    :known => [x2, x3],
    :correct => [a, b + c, x1, x2, x3],
))


@testset "Identifiable functions with known generic initial conditions" begin
    for case in cases
        ode = case[:ode]
        known = case[:known]
        result = find_identifiable_functions_kic(ode, known)
        correct = [f // one(parent(ode)) for f in case[:correct]]
        @test Set(result) == Set(correct)
    end
end
