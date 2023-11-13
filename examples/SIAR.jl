# SIAR model
# L. Gallo, M.Frasca, V. Latora, G. Russo
# Lack of practical identifiability may hamper reliable predictions in COVID-19 epidemic models
# Science Advances, DOI: 10.1126/sciadv.abg5234
# Equation (21)
using StructuralIdentifiability

ode = @ODEmodel(
    S'(t) =
        -S(t) * (betaIa * Ia(t) + betaIs * Is(t) + betaH * H(t) + betaT * T(t)) * Ninv,
    E'(t) =
        S(t) * (betaIa * Ia(t) + betaIs * Is(t) + betaH * H(t) + betaT * T(t)) * Ninv -
        (alphaEIa + alphaEIs) * E(t),
    Ia'(t) = alphaEIa * E(t) - (alphaIaIs + alphaIaRu) * Ia(t) - xi * Ia(t),
    Is'(t) =
        alphaEIs * E(t) + alphaIaIs * Ia(t) -
        (alphaIsH + alphaIsT + alphaIsRu + alphaIsD) * Is(t),
    H'(t) = alphaIsH * Is(t) + xi * Ia(t) - (alphaHT + alphaHRd) * H(t),
    T'(t) = alphaIsT * Is(t) + alphaHT * H(t) - (alphaTRd + alphaTD) * T(t),
    Rd'(t) = alphaHRd * H(t) + alphaTRd * T(t),
    D'(t) = alphaIsD * Is(t) + alphaTD * T(t),
    y1(t) = H(t),
    y2(t) = T(t),
    y3(t) = Rd(t),
    y4(t) = D(t),
    y5(t) = Ninv
)

println(find_identifiable_functions(ode, funcs_to_check = ode.parameters))
