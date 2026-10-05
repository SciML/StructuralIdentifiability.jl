using StructuralIdentifiability
using StructuralIdentifiability: cmp_states, cmp_lie_states

ode = @ODEmodel(
    A'(t) = a * A(t),
    B'(t) = b * B(t),
    C'(t) = c * C(t),
    y(t) = A(t) + B(t) + C(t)
)

println(observation_field(ode, cmp=cmp_states(ode)))
println(observation_field(ode, cmp=cmp_lie_states(ode)))
