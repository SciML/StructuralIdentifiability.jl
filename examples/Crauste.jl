using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(
    N'(t) = -N(t) * mu_N - N(t) * P(t) * delta_NE,
    E'(t) =
        N(t) * P(t) * delta_NE - E(t)^2 * mu_EE - E(t) * delta_EL + E(t) * P(t) * rho_E,
    S'(t) = S(t) * delta_EL - S(t) * delta_LM - S(t)^2 * mu_LL - E(t) * S(t) * mu_LE,
    M'(t) = S(t) * delta_LM - mu_M * M(t),
    P'(t) = P(t)^2 * rho_P - P(t) * mu_P - E(t) * P(t) * mu_PE - S(t) * P(t) * mu_PL,
    y1(t) = N(t),
    y2(t) = E(t) + S(t),
    y3(t) = M(t)
)

@time println(assess_identifiability(ode))
