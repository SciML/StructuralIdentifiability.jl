# https://doi.org/10.1016/j.jtbi.2014.06.020
# Output y1 = x chosen based on discussion in Section 4.2

# Eq. 2
ode = @ODEmodel(
    x'(t) = rho1*x(t)*(1 - x(t)) - alpha13*x(t)*z(t),
    y'(t) = rho2*y(t)*z(t)/(1 + z(t)) - alpha23*y(t)*z(t) - delta2*y(t),
    z'(t) = rho3*z(t)*(1 - z(t)) - alpha31*z(t)*x(t) - alpha32*z(t)*y(t),
    y1(t) = x(t)
)

# OK!
assess_identifiability(ode)

# Eq. 7
ode = @ODEmodel(
    x'(t) = rho1*x(t)*(1 - x(t)) - alpha13*x(t)*z(t),
    y'(t) = rho2*y(t)*z(t)/(1 + z(t)) - alpha23*y(t)*z(t) - delta2*y(t) + alpha24*y(t)*w(t),
    z'(t) = rho3*z(t)*(1 - z(t)) - alpha31*z(t)*x(t) - alpha32*z(t)*y(t) + alpha34*z(t)*w(t)/(1 + w(t)),
    w'(t) = rho4*w(t)*z(t)/(1 + z(t)) - delta4*w(t),
    y1(t) = x(t)
)

### Seems hard
assess_identifiability(ode)

### Seems hard
include(joinpath(@__DIR__, "shared.jl"))
gens, simple_gens = simplify_lie(ode, 7);
println(simple_gens)
