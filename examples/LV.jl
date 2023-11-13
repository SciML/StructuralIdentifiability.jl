# Lotka-Volterra model
using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = a * x1(t) - b * x1(t) * x2(t) + u(t),
    x2'(t) = -c * x2(t) + d * x1(t) * x2(t),
    y(t) = x1(t)
)

println(assess_identifiability(ode))

# Not everything is identifiable, so we may wonder what are the identifiable functions

println(find_identifiable_functions(ode, with_states = true))
