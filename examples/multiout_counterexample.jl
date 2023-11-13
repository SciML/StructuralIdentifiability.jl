# an artificial example illustrating the necessity of an extra projection in the multi-output case
using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = (1 + x1(t)^2) // 2,
    x2'(t) = (1 - x1(t)^2) // (1 + x1(t)^2),
    y1(t) = 2 * x1(t) // (b * (1 + x1(t)^2)),
    y2(t) = x2(t)
)

@time println(assess_identifiability(ode))
