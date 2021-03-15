using Logging

include("../src/StructuralIdentifiability.jl")
using .StructuralIdentifiability

# model 41, normalized (all state variables are divided by N)
logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(
    S'(t) = -b * S(t) * (I(t) + J(t) + q * A(t)),
    E'(t) = b * S(t) * (I(t) + J(t) + q * A(t)) - k * E(t),
    A'(t) = k * (1 - r) * E(t) - g1 * A(t),
    I'(t) = k * r * E(t) - (alpha + g1) * I(t),
    J'(t) = alpha * I(t) - g2 * J(t),
    C'(t) = alpha * I(t),
    y(t) = C(t)
)

@time println(assess_global_identifiability(ode))
