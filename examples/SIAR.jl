# SIAR model
# https://arxiv.org/pdf/2012.00443.pdf
using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

ode = @ODEmodel(function S'(t)
                    -S(t) *
                    (betaIa * Ia(t) + betaIs * Is(t) + betaH * H(t) + betaT * T(t)) * Ninv
                end,
                function E'(t)
                    S(t) * (betaIa * Ia(t) + betaIs * Is(t) + betaH * H(t) + betaT * T(t)) *
                    Ninv - (alphaEIa + alphaEIs) * E(t)
                end,
                Ia'(t)=alphaEIa * E(t) - (alphaIaIs + alphaIaRu) * Ia(t) - xi * Ia(t),
                function Is'(t)
                    alphaEIs * E(t) + alphaIaIs * Ia(t) -
                    (alphaIsH + alphaIsT + alphaIsRu + alphaIsD) * Is(t)
                end,
                H'(t)=alphaIsH * Is(t) + xi * Ia(t) - (alphaHT + alphaHRd) * H(t),
                T'(t)=alphaIsT * Is(t) + alphaHT * H(t) - (alphaTRd + alphaTD) * T(t),
                Rd'(t)=alphaHRd * H(t) + alphaTRd * T(t),
                D'(t)=alphaIsD * Is(t) + alphaTD * T(t),
                y1(t)=H(t),
                y2(t)=T(t),
                y3(t)=Rd(t),
                y4(t)=D(t),
                y5(t)=Ninv)

@time println(assess_identifiability(ode))
