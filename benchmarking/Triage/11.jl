# https://doi.org/10.1016/S0025-5564(97)10004-9

using StructuralIdentifiability

# eq. 16-18
ode = @ODEmodel(
    x'(t) = p1*x(t)/(p2 + p3*x(t)) + u(t),
    y(t) = c*x(t)
)

find_identifiable_functions(ode, with_states=true)

# eq. 44-47
ode2 = @ODEmodel(
    x'(t) = mum*s(t)*x(t)/(Ks + s(t)) - Kd*x(t) + b1*u(t),
    s'(t) = -mum*s(t)*x(t)/(Y*(Ks + s(t))) + b2*u(t),
    y(t) = x(t)
)

find_identifiable_functions(ode, with_states=true)
