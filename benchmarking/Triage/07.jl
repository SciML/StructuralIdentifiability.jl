# https://doi.org/10.1103/PhysRevE.89.050902

using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = -x1(t)^3 + p1*x1(t)^2 + x2(t) - x3(t),
    x2'(t) = 1 - p2*x1(t)^2 - x2(t),
    x3'(t) = p3*(x1(t) + p4*(p5 - x3(t))),
    y(t) = x1(t)
)

id = find_identifiable_functions(ode; with_states=true)
