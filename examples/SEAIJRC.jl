# K. Roosa and G. Chowell. 
# Assessing parameter identifiability in compartmental dynamic models using a computational approach: application to infectious disease transmission models
# https://doi.org/10.1186/s12976-018-0097-6
# Model 2
using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(S'(t) = -b * S(t) * (I(t) + J(t) + q * A(t)) * Ninv,
                E'(t) = b * S(t) * (I(t) + J(t) + q * A(t)) * Ninv - k * E(t),
                A'(t) = k * (1 - r) * E(t) - g1 * A(t),
                I'(t) = k * r * E(t) - (alpha + g1) * I(t),
                J'(t) = alpha * I(t) - g2 * J(t),
                C'(t) = alpha * I(t),
                y(t) = C(t),
                y2(t) = Ninv)

@time println(assess_global_identifiability(ode))
