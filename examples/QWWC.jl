using StructuralIdentifiability

ode = @ODEmodel(
    x'(t) = a * (y(t) - x(t)) + y(t) * z(t),
    y'(t) = b * (x(t) + y(t)) - x(t) * z(t),
    z'(t) = -c * z(t) - d * w(t) + x(t) * y(t),
    w'(t) = e * z(t) - f * w(t) + x(t) * y(t),
    g(t) = x(t)
)

println(assess_identifiability(ode))
