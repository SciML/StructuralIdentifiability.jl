using Revise, StructuralIdentifiability, ParamPunPam, RationalFunctionFields
import StructuralIdentifiability: parent_ring_change, lie_derivative
import RationalFunctionFields: fields_equal

cmp = (f, g) -> begin
    f = parent_ring_change(f, parent(ode))
    g = parent_ring_change(g, parent(ode))
    df = lie_derivative(f, ode)
    dg = lie_derivative(g, ode)
    default_cmp = RationalFunctionFields.rational_function_cmp
    if iszero(df) && iszero(dg)
        default_cmp(f, g)
    else
        default_cmp(df, dg)
    end    
end

ode = @ODEmodel(
            x1'(t) =
                (p1 * x4(t)) - (p3 * x1(t)) -
                p4 * ((x1(t)^2 / (p5 + x1(t))) * (1 + (p6 * u1(t) / (p7 + u1(t))))),
            x2'(t) =
                p8 - (p9 * x2(t)) -
                p10 * (
                (x1(t) * x2(t) / (p11 + x2(t))) * (1 + (p12 * u1(t) / (p13 + u1(t))))
            ),
            x3'(t) =
                p14 - (p15 * x3(t)) -
                p16 * x1(t) * x3(t) * (1 - p18 * u1(t)) / (p17 + x3(t)),
            x4'(t) =
                p20 - p21 * (1 - p24) * (1 - p25) / ((p22^4) + 1) - (p20 * x4(t)) +
                (p21 * (x3(t)^4)) +
                (1 + p23 * u1(t)) * (1 - p24 * x1(t)) * (1 - p25 * x2(t)) /
                (p22^4 + x3(t)^4),
            y1(t) = x1(t),
            y2(t) = x2(t),
            y3(t) = x3(t),
            y4(t) = x4(t)
        )

id1 = find_identifiable_functions(ode, with_states=true, cmp = cmp)

id2 = find_identifiable_functions(ode, with_states=true, cmp = RationalFunctionFields.rational_function_cmp)

@assert fields_equal(RationalFunctionField(id1), RationalFunctionField(id2), 0.99)
