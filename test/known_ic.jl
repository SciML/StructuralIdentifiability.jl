cases = []

ode = @ODEmodel(
    x1'(t) = a * x1(t) - b * x1(t) * x2(t),
    x2'(t) = -c * x2(t) + d * x1(t) * x2(t),
    y(t) = x1(t)
)

push!(cases, Dict(
    :ode => ode,
    :known => [x2],
    :to_check => [],
    :correct_funcs => [a, b, c, d, x1, x2],
    :correct_ident => OrderedDict(x => :globally for x in [x1 , x2, a, b, c, d]),
))

ode = @ODEmodel(
    x1'(t) = a  + x2(t) + x3(t),
    x2'(t) = b^2 + c,
    x3'(t) = -c,
    y(t) = x1(t)
)

push!(cases, Dict(
    :ode => ode,
    :known => [x2, x3],
    :to_check => [],
    :correct_funcs => [a, b^2, x1, x2, x3],
    :correct_ident => OrderedDict(x1 => :globally, x2 => :globally, x3 => :globally, a => :globally, b => :locally, c => :nonidentifiable),
))

push!(cases, Dict(
    :ode => ode,
    :known => [x2, x3],
    :to_check => [b^2, x2 * c],
    :correct_funcs => [a, b^2, x1, x2, x3],
    :correct_ident => OrderedDict(b^2 => :globally, x2 * c => :nonidentifiable),
))

@testset "Identifiable functions with known generic initial conditions" begin
    for case in cases
        ode = case[:ode]
        known = case[:known]

        result_funcs = find_identifiable_functions_kic(ode, known)
        correct_funcs = [f // one(parent(ode)) for f in case[:correct_funcs]]
        @test Set(result_funcs) == Set(correct_funcs)

        result_ident = assess_identifiability_kic(ode, known, funcs_to_check = case[:to_check])
        @test case[:correct_ident] == result_ident
    end
end
