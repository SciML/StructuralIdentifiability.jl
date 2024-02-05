cases = []

ode = @ODEmodel(
    x1'(t) = a * x1(t) - b * x1(t) * x2(t),
    x2'(t) = -c * x2(t) + d * x1(t) * x2(t),
    y(t) = x1(t)
)

push!(
    cases,
    Dict(
        :ode => ode,
        :known => [x2],
        :to_check => [],
        :correct_funcs => [a, b, c, d, x1, x2],
        :correct_ident => OrderedDict(x => :globally for x in [x1, x2, a, b, c, d]),
    ),
)

ode = @ODEmodel(x1'(t) = a + x2(t) + x3(t), x2'(t) = b^2 + c, x3'(t) = -c, y(t) = x1(t))

push!(
    cases,
    Dict(
        :ode => ode,
        :known => [x2, x3],
        :to_check => [],
        :correct_funcs => [a, b^2, x1, x2, x3],
        :correct_ident => OrderedDict(
            x1 => :globally,
            x2 => :globally,
            x3 => :globally,
            a => :globally,
            b => :locally,
            c => :nonidentifiable,
        ),
    ),
)

push!(
    cases,
    Dict(
        :ode => ode,
        :known => [x2, x3],
        :to_check => [b^2, x2 * c],
        :correct_funcs => [a, b^2, x1, x2, x3],
        :correct_ident => OrderedDict(b^2 => :globally, x2 * c => :nonidentifiable),
    ),
)

ode = @ODEmodel(x1'(t) = a * x1(t), x2'(t) = x1(t) + 1 / x1(t), y(t) = x2(t))

push!(
    cases,
    Dict(
        :ode => ode,
        :known => [x1],
        :to_check => [],
        :correct_funcs => [a, x1, x2],
        :correct_ident => OrderedDict(x1 => :globally, x2 => :globally, a => :globally),
    ),
)

ode = @ODEmodel(
    x1'(t) = -b * x1(t) + 1 / (c + x4(t)),
    x2'(t) = alpha * x1(t) - beta * x2(t),
    x3'(t) = gama * x2(t) - delta * x3(t),
    x4'(t) = sigma * x4(t) * (gama * x2(t) - delta * x3(t)) / x3(t),
    y(t) = x1(t)
)

push!(
    cases,
    Dict(
        :ode => ode,
        :known => [x2, x3],
        :to_check => [alpha, alpha * gama],
        :correct_funcs => [
            sigma,
            c,
            b,
            x4,
            x3,
            x2,
            x1,
            beta * delta,
            alpha * gama,
            beta + delta,
            -delta * x3 + gama * x2,
        ],
        :correct_ident => OrderedDict(alpha => :locally, alpha * gama => :globally),
    ),
)

@testset "Identifiable functions with known generic initial conditions" begin
    for case in cases
        ode = case[:ode]
        known = case[:known]

        result_funcs = find_identifiable_functions(ode, known_ic = known)
        correct_funcs =
            replace_with_ic(ode, [f // one(parent(ode)) for f in case[:correct_funcs]])
        @test Set(result_funcs) == Set(correct_funcs)

        result_ident =
            assess_identifiability(ode, known_ic = known, funcs_to_check = case[:to_check])
        @test OrderedDict(
            replace_with_ic(ode, [k])[1] => v for (k, v) in case[:correct_ident]
        ) == result_ident
    end
end
