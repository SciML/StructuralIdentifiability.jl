using Logging

include("../src/StructuralIdentifiability.jl")
using .StructuralIdentifiability

# SIAR model
# https://arxiv.org/pdf/2012.00443.pdf
logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(
    S'(t) = -S(t) * (betaIa * Ia(t) + betaIs * Is(t) + betaH * H(t) + betaT * T(t)) * Ninv(t),
    E'(t) = S(t) * (betaIa * Ia(t) + betaIs * Is(t) + betaH * H(t) + betaT * T(t)) * Ninv(t) - (alphaEIa + alphaEIs) * E(t),
    Ia'(t) = alphaEIa * E(t) - (alphaIaIs + alphaIaRu) * Ia(t) - xi * Ia(t),
    Is'(t) = alphaEIs * E(t) + alphaIaIs * Ia(t) - (alphaIsH + alphaIsT + alphaIsRu + alphaIsD) * Is(t),
    H'(t) = alphaIsH * Is(t) + xi * Ia(t) - (alphaHT + alphaHRd) * H(t),
    T'(t) = alphaIsT * Is(t) + alphaHT * H(t) - (alphaTRd + alphaTD) * T(t),
    Rd'(t) = alphaHRd * H(t) + alphaTRd * T(t),
    D'(t) = alphaIsD * Is(t) + alphaTD * T(t),
    Ninv'(t) = 0,
    y1(t) = H(t),
    y2(t) = T(t),
    y3(t) = Rd(t),
    y4(t) = D(t),
    y5(t) = Ninv(t)
)

@time println(assess_identifiability(ode))

#@time ifunc = find_identifiable_functions(ode)

#println(ifunc)
