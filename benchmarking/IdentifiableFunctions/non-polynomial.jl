using StructuralIdentifiability

# 1.
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
@info "" funcs
#=
Identifiabile functions:
│     I
│     z
│     d
│     w*E
│     w*S
│     w*U + w*I
└     beta//(w*N)
=#

# 2.
SLIQR = StructuralIdentifiability.@ODEmodel(
    S'(t) = -b * In(t) * S(t) * Ninv - u(t) * S(t) * Ninv,
    L'(t) = b * In(t) * S(t) * Ninv - a * L(t),
    In'(t) = a * L(t) - g * In(t) + s * Q(t),
    Q'(t) = (1 - e) * g * In(t) - s * Q(t),
    y(t) = In(t) * Ninv
)
funcs1 = StructuralIdentifiability.find_identifiable_functions(
    SLIQR,
    with_states = true,
    strategy = (:normalforms, 2),
    seed = 42,
)

StructuralIdentifiability.field_contains(
    StructuralIdentifiability.RationalFunctionField(funcs1),
    [(s * Q - Q * a) // one(s)],
    0.99,
)

funcs2 = StructuralIdentifiability.find_identifiable_functions(
    SLIQR,
    with_states = true,
    strategy = (:hybrid, 3),
    seed = 42,
)
@info "" funcs1
@info "" funcs2

# nf(y6 y10 - C1) = 932277181*y10 + 573449642
# nf(y4 y10 - C2) = 932277181*y10 + 573449642

R, (s, a, Q) = QQ["s", "a", "Q"]
f = [s, Q^2, s * Q - a * Q]

R, (y1, y2, y3) = QQ["y1", "y2", "y3"]
StructuralIdentifiability.relations_over_ff([y1 - 1, y1 + 10, y1 + 10], [s, a, Q])

rff = StructuralIdentifiability.RationalFunctionField(f)
relations = StructuralIdentifiability.monomial_generators_up_to_degree(rff, 2)

StructuralIdentifiability.fields_equal(
    StructuralIdentifiability.RationalFunctionField(funcs1),
    StructuralIdentifiability.RationalFunctionField(funcs2),
    0.99,
)

@info "" funcs
#=
│     In
│     s
│     Ninv
│     b
│     S*a
│     e*g*a
│     g + a
│     s*Q - Q*a
│     e*s*g - s*g + g*a
└     (In + Q + L)//S
=#
