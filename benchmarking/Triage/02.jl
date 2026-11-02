# https://doi.org/10.1016/j.cmpb.2013.06.013
# eq. 43-48

using StructuralIdentifiability

ode = @ODEmodel(
        x1'(t) = k4*x4(t) - r4*x1(t) - k1*x1(t)*(T0 - x2(t)) + r1*x2(t),

        x2'(t) = k1*x1(t)*(T0 - x2(t)) - (r1 + k2)*x2(t),

        x3'(t) = r3*x4(t) - k3*x3(t) + k2*x2(t),

        x4'(t) = r4*x1(t) - (k4 + r3)*x4(t) + k3*x3(t),

        y(t) = x2(t) + x3(t) + x4(t)
    )

find_identifiable_functions(ode; with_states = true)
