using ModelingToolkit
using StructuralIdentifiability

@parameters α β
@variables t S(t) I(t) R(t) y(t)
D = Difference(t; dt = 1.0)

eqs = [D(S) ~ S - β * S * I, D(I) ~ I + β * S * I - α * I, D(R) ~ R + α * I]
@named sir = DiscreteSystem(eqs)

println(assess_local_identifiability(sir; measured_quantities = [y ~ I]))
