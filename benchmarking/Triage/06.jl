# https://doi.org/10.1098/rsta.2011.0544
# Supplementary Notes

using StructuralIdentifiability, Nemo

ode = @ODEmodel(
    x1'(t) = -kon*x1(t)*x2(t) + kon*kD*x3(t) + kex*x4(t),
    x2'(t) = -kon*x1(t)*x2(t) + kon*kD*x3(t) + kt*Bmax - kt*x2(t) + kex*x4(t),
    x3'(t) = kon*x1(t)*x2(t) - kon*kD*x3(t) - ke*x3(t),
    x4'(t) = ke*x3(t) - (kex + kdi + kde)*x4(t),
    x5'(t) = kdi*x4(t),
    x6'(t) = kde*x4(t),

    y1(t) = scale*(x1(t) + x6(t)),
    y2(t) = scale*x3(t),
    y3(t) = scale*(x4(t) + x5(t))
)

id = find_identifiable_functions(ode; with_states=true)
