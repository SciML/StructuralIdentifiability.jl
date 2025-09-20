# Test that the following pipeline works:
# ODE -> IO equations -> raw identifiable functions -> RFF -> MQS -> (maybe specialize) -> ideal generators
import Groebner

@testset "Raw Generators of RFF" begin
    # SIWR
    ode = StructuralIdentifiability.@ODEmodel(
        S'(t) = mu - bi * S(t) * I(t) - bw * S(t) * W(t) - mu * S(t) + a * R(t),
        I'(t) = bw * S(t) * W(t) + bi * S(t) * I(t) - (gam + mu) * I(t),
        W'(t) = xi * (I(t) - W(t)),
        R'(t) = gam * I(t) - (mu + a) * R(t),
        y(t) = k * I(t)
    )

    io_eqs = StructuralIdentifiability.find_ioequations(ode)
    id_funcs, bring = StructuralIdentifiability.extract_identifiable_functions_raw(
        io_eqs,
        ode,
        empty(ode.parameters),
        true,
    )

    param_ring, _ = polynomial_ring(
        base_ring(bring),
        map(string, ode.parameters),
        internal_ordering = Nemo.internal_ordering(bring),
    )

    id_funcs_no_states = map(
        polys -> map(
            poly -> StructuralIdentifiability.parent_ring_change(poly, param_ring),
            polys,
        ),
        id_funcs[:no_states],
    )

    rff = RationalFunctionField(id_funcs_no_states)

    # Part 1: mod p and specialized
    p = Nemo.Native.GF(2^62 + 135)
    StructuralIdentifiability.ParamPunPam.reduce_mod_p!(rff.mqs, p)
    point = rand(
        p,
        length(Nemo.gens(StructuralIdentifiability.ParamPunPam.parent_params(rff.mqs))),
    )
    eqs = StructuralIdentifiability.ParamPunPam.specialize_mod_p(rff.mqs, point)
    gb = Groebner.groebner(eqs, ordering = Groebner.DegRevLex())
    # GB is linear
    @test length(gb) == length(gens(parent(eqs[1])))
    expected = 9202476
    str = join(map(string, eqs), ",")
    @info "" length(str)
    @test abs(length(str) - expected) / expected * 100 < 5

    # Part 2: over Q
    eqs = fractionfree_generators_raw(rff.mqs)[1]
    expected = 21486079
    str = join(map(string, eqs), ",")
    @info "" length(str)
    @test abs(length(str) - expected) / expected * 100 < 5
end
