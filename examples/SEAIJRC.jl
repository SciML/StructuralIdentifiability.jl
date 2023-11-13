# K. Roosa and G. Chowell. 
# Assessing parameter identifiability in compartmental dynamic models using a computational approach: application to infectious disease transmission models
# https://doi.org/10.1186/s12976-018-0097-6
# Model 2
using StructuralIdentifiability

ode = @ODEmodel(
    S'(t) = -b * S(t) * (I(t) + J(t) + q * A(t)) * Ninv,
    E'(t) = b * S(t) * (I(t) + J(t) + q * A(t)) * Ninv - k * E(t),
    A'(t) = k * (1 - r) * E(t) - g1 * A(t),
    I'(t) = k * r * E(t) - (alpha + g1) * I(t),
    J'(t) = alpha * I(t) - g2 * J(t),
    C'(t) = alpha * I(t),
    y(t) = C(t),
    y2(t) = Ninv
)

println(assess_identifiability(ode))

# Not everything is identifiable, so we may wonder what are the identifiable functions

println(find_identifiable_functions(ode, with_states = true))
