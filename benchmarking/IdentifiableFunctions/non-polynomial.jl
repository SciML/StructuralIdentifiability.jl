using StructuralIdentifiability

# 1.
# Non-polynomial
SEUIR = StructuralIdentifiability.@ODEmodel(
    S'(t) = -beta * (U(t) + I(t)) * (S(t) / N),
    E'(t) = beta * (U(t) + I(t)) * (S(t) / N) - E(t) * z,
    U'(t) = (z - w) * E(t) - U(t) * d,
    I'(t) = w * E(t) - I(t) * d,
    R'(t) = (U(t) + I(t)) * d,
    y1(t) = I(t)
)
funcs = StructuralIdentifiability.find_identifiable_functions(
    SEUIR,
    with_states = true,
    strategy = (:hybrid, 3),
    seed = 42,
)

# 2.
# Goodwin oscillator
# ?
Goodwin = StructuralIdentifiability.@ODEmodel(
    x1'(t) = -b * x1(t) + 1 / (c + x4(t)),
    x2'(t) = alpha * x1(t) - beta * x2(t),
    x3'(t) = gama * x2(t) - delta * x3(t),
    x4'(t) = sigma * x4(t) * (gama * x2(t) - delta * x3(t)) / x3(t),
    y(t) = x1(t)
)
funcs = StructuralIdentifiability.find_identifiable_functions(
    Goodwin,
    with_states = true,
    strategy = (:hybrid, 3),
    seed = 42,
)

# 3.
# KD1999
# Non-polynomial
KD1999 = StructuralIdentifiability.@ODEmodel(
    Ca'(t) = u1(t) * (Ca0 - Ca(t)) / V - k0 * Arr * Ca(t),
    Cb'(t) = -u1(t) * Cb(t) / V + k0 * Arr(t) * Ca(t),
    T'(t) =
        u1(t) * (Ta - T(t)) / V -
        (k0 * Arr(t) * Ca(t) * DH + UA * (Tj(t) - T(t)) / V) / (ro * cp),
    Tj'(t) = u2(t) * (Th - Tj(t)) / Vh - UA / (roh * cph) * (Tj(t) - T(t)) / Vh,
    Arr'(t) =
        E * Arr(t) / (R * T(t)^2) * (
            u1(t) * (Ta - T(t)) / V -
            (k0 * Arr(t) * Ca(t) * DH + UA * (Tj(t) - T(t)) / V) / (ro * cp)
        ),
    y1(t) = Cb(t),
    y2(t) = T(t)
)
funcs = StructuralIdentifiability.find_identifiable_functions(
    KD1999,
    with_states = true,
    strategy = (:normalforms, 3),
    seed = 42,
)

# 4.
# LLW1987_io
# Non-polynomial
LLW1987_io = StructuralIdentifiability.@ODEmodel(
    x1'(t) = -p1 * x1(t) + p2 * u(t),
    x2'(t) = -p3 * x2(t) + p4 * u(t),
    x3'(t) = -(p1 + p3) * x3(t) + (p4 * x1(t) + p2 * x2(t)) * u(t),
    y1(t) = x3(t)
)
funcs = StructuralIdentifiability.find_identifiable_functions(
    LLW1987_io,
    with_states = true,
    strategy = (:normalforms, 3),
    seed = 42,
)
