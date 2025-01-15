# For each ODE system we check the equality (in terms of fields of rational
# functions) of the true set of identifiable functions and the obtained
# simplified set
test_cases = []

###

ode = StructuralIdentifiability.@ODEmodel(x'(t) = a * x(t) + u(t), y(t) = x(t))
ident_funcs = [a]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

ode = StructuralIdentifiability.@ODEmodel(x1'(t) = a, x2'(t) = -a, y(t) = x1 + x2)
ident_funcs = []
push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# Parameter a is not identifiable, and neither are any combinations thereof.
ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = x2(t) - a,
    x2'(t) = x1(t) + a,
    y(t) = x1(t) + x2(t)
)
ident_funcs =
    Vector{StructuralIdentifiability.AbstractAlgebra.Generic.FracFieldElem{typeof(x1)}}()

push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# Example 2 from 
# "On Global Identifiability for Arbitrary Model Parametrizations",
# DOI: 10.1016/0005-1098(94)90029-9
ode = StructuralIdentifiability.@ODEmodel(x1'(t) = Θ * x2(t)^2, x2'(t) = u(t), y(t) = x1(t))
ident_funcs = [Θ]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# Example 4 from 
# "On Global Identifiability for Arbitrary Model Parametrizations",
# DOI: 10.1016/0005-1098(94)90029-9
ode = StructuralIdentifiability.@ODEmodel(
    x'(t) = (-V_m * x(t)) / (k_m + x(t)) + k01 * x(t),
    y(t) = c * x(t)
)
ident_funcs = [k01, c * k_m, V_m * c]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# Parameters b and c enter the io equations only as the product b * c.
ode = StructuralIdentifiability.@ODEmodel(x'(t) = a * x(t) + b * u(t), y(t) = c * x(t))
ident_funcs = [b * c, a]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# 2-compartiment model
ode = StructuralIdentifiability.@ODEmodel(
    x0'(t) = -(a01 + a21) * x0(t) + a12 * x1(t),
    x1'(t) = a21 * x0(t) - a12 * x1(t),
    y(t) = x0(t)
)
ident_funcs = [(a01 * a12), (a01 + a12 + a21)]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# TODO: uncomment when identifiability can handle models with no states 
# ode = StructuralIdentifiability.@ODEmodel(
#     y(t) = a*u(t)
# )
# ident_funcs = [(a01 * a12), (a01 + a12 + a21)]
# push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# Example 2 from
# "On Structural Identifiability",
# DOI: https://doi.org/10.1016/0025-5564(70)90132-X
#
# More or less the same 2-compartmental model as the one given above
ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = -(k1 + k2) * x1(t) + k3 * x2(t) + u(t),
    x2'(t) = k2 * x1(t) - (k3 + k4) * x2(t),
    y(t) = x1(t)
)
ident_funcs = [k1 + k2, k3 + k4, k2 * k3]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# Diagonal with simple spectrum and observable states
ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = λ1 * x1(t) + β1 * u1(t),
    x2'(t) = λ2 * x2(t) + β2 * u2(t),
    x3'(t) = λ3 * x3(t) + β3 * u3(t),
    y(t) = x1(t) + x2(t) + x3(t)
)
ident_funcs = [λ1, λ2, λ3, β1, β2, β3]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# 3 compartments:
#   x1 <--> x2 <--> x3
# If we observe x1 and control x1, then all parameters are identifiable
ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = -a1 * x1(t) + b1 * x2(t) + u(t),
    x2'(t) = -(a2 + b1) * x2(t) + a1 * x1(t) + b2 * x3(t),
    x3'(t) = -b2 * x3(t) + a2 * x2(t),
    y(t) = x1(t)
)
ident_funcs = [a1, a2, b1, b2]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# Example 3 from
# "On Structural Identifiability",
# DOI: https://doi.org/10.1016/0025-5564(70)90132-X
# 
# 3 compartments:
#   x1 <--> x2 <--> x3
# If we observe x1 and control x3, then only some functions of parameters
# are identifiable
ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = -a1 * x1(t) + b1 * x2(t),
    x2'(t) = -(a2 + b1) * x2(t) + a1 * x1(t) + b2 * x3(t),
    x3'(t) = -b2 * x3(t) + a2 * x2(t) + u(t),
    y(t) = x1(t)
)
ident_funcs = [b1 * b2, a1 + a2 + b1 + b2, a1 * a2 + a1 * b2]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# Example 3 from
# "On the identifiability and distinguishability of nonlinear parametric
# models",
# DOI: https://doi.org/10.1016/0378-4754(95)00123-9
ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = p1 * x1^2 + p2 * x1 * x2,
    x2'(t) = p3 * x1^2 + p4 * x1 * x2,
    y(t) = x1
)
ident_funcs = [p1 + p4, -p2 * p3 + p4 * p1]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# Goowdin oscillator
ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = -b * x1(t) + 1 / (c + x4(t)),
    x2'(t) = alpha * x1(t) - beta * x2(t),
    x3'(t) = gama * x2(t) - delta * x3(t),
    x4'(t) = sigma * x4(t) * (gama * x2(t) - delta * x3(t)) / x3(t),
    y(t) = x1(t)
)
ident_funcs = [sigma, delta + beta, c, b, delta * beta]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# SIRS forced
ode = StructuralIdentifiability.@ODEmodel(
    s'(t) = mu - mu * s(t) - b0 * (1 + b1 * x1(t)) * i(t) * s(t) + g * r(t),
    i'(t) = b0 * (1 + b1 * x1(t)) * i(t) * s(t) - (nu + mu) * i(t),
    r'(t) = nu * i(t) - (mu + g) * r(t),
    x1'(t) = -M * x2(t),
    x2'(t) = M * x1(t),
    y1(t) = i(t),
    y2(t) = r(t)
)
ident_funcs = [g, mu, b0, nu, M^2]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# SEIR_1_io
ode = StructuralIdentifiability.@ODEmodel(
    S'(t) = -beta * S(t) * I(t),
    E'(t) = beta * S(t) * I(t) - v * E(t),
    I'(t) = v * E(t) - psi * I(t) - (1 - psi) * gamma * I(t),
    R'(t) = gamma * Q(t) + (1 - psi) * gamma * I(t),
    Q'(t) = -gamma * Q(t) + psi * I(t),
    y1(t) = Q(t)
)
ident_funcs = [gamma * psi - psi * v, beta // psi, gamma, psi * v - psi - v]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# Bilirubin2_io.
# Regression test: failed before, as the total degrees were being estimated
# incorrectly in the interpolation
ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) =
        -(k21 + k31 + k41 + k01) * x1(t) + k12 * x2(t) + k13 * x3(t) + k14 * x4(t) + u(t),
    x2'(t) = k21 * x1(t) - k12 * x2(t),
    x3'(t) = k31 * x1(t) - k13 * x3(t),
    x4'(t) = k41 * x1(t) - k14 * x4(t),
    y1(t) = x1(t)
)
ident_funcs = [
    k01 // one(k01),
    k12 * k13 * k14 // one(k01),
    k31 * k21 * k41 // one(k01),
    k12 + k13 + k14 // one(k01),
    k31 + k21 + k41 // one(k01),
    k12 * k13 + k12 * k14 + k13 * k14 // one(k01),
    k31 * k21 + k31 * k41 + k21 * k41 // one(k01),
    k31 * k12 - 2 * k31 * k13 + k31 * k14 - 2 * k21 * k12 +
    k21 * k13 +
    k21 * k14 +
    k12 * k41 +
    k13 * k41 - 2 * k14 * k41 // one(k01),
]
# Too slow with hybrid strategy :(
# push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# Biohydrogenation_io
ode = StructuralIdentifiability.@ODEmodel(
    x5'(t) =
        (
            k5 * k8 * x4(t) + k5 * x6(t) * x4(t) + k5 * x5(t) * x4(t) - k6 * x5(t) * k7 - x5(t) * k7 * x4(t)
        ) // (
            k8 * k6 + k8 * x4(t) + k6 * x6(t) + k6 * x5(t) + x6(t) * x4(t) + x5(t) * x4(t)
        ),
    x7'(t) = (k9 * k10 * x6(t) - k9 * x6(t)^2) // k10,
    x4'(t) = (-k5 * x4(t)) // (k6 + x4(t)),
    x6'(t) =
        (
            -k8 * k9 * k10 * x6(t) + k8 * k9 * x6(t)^2 - k9 * k10 * x6(t)^2 -
            k9 * k10 * x6(t) * x5(t) +
            k9 * x6(t)^3 +
            k9 * x6(t)^2 * x5(t) +
            k10 * x5(t) * k7
        ) // (k8 * k10 + k10 * x6(t) + k10 * x5(t)),
    y1(t) = x4(t),
    y2(t) = x5(t)
)
ident_funcs = [k7, k5, k6, k10 * k9, k9^2, k10 + 2 * k8]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# SLIQR
ode = StructuralIdentifiability.@ODEmodel(
    S'(t) = -b * In(t) * S(t) * Ninv - S(t) * Ninv * u(t),
    In'(t) = -In(t) * g + s * Q(t) + a * L(t),
    L'(t) = b * In(t) * S(t) * Ninv - a * L(t),
    Q'(t) = -e * In(t) * g + In(t) * g - s * Q(t),
    y(t) = In(t) * Ninv
)
ident_funcs = [
    (a * e) // (a + e * s - s),
    b,
    a + g,
    (
        a^2 * e * s + a^2 * g + 3 * a * e * g * s - a * e * s^2 - 2 * a * g * s +
        e^2 * g * s^2 - 2 * e * g * s^2 + g * s^2
    ) // (a + e * s - s),
    s,
    Ninv,
]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# St.
# Regression test:
# Failed before, as the degrees of Groebner basis were too large
ode = StructuralIdentifiability.@ODEmodel(
    S'(t) = -e * S(t) - S(t) * d * W(t) + S(t) * r - S(t) * a * W(t) + R(t) * g,
    R'(t) = e * S(t) + rR * R(t) + S(t) * a * W(t) - dr * R(t) * W(t) - R(t) * g,
    W'(t) = T * Dd - W(t) * Dd,
    y1(t) = S(t) + R(t),
    y2(t) = T
)
ident_funcs = [
    T,
    Dd,
    e - rR + dr * T + d * T + g - r + a * T,
    (d * rR - dr * r) // (d - dr),
    (dr^2 + d^2 + 2 * d * a + a^2) // (dr * d + dr * a),
    (e * dr - e * d + rR * a + dr * g - d * g - r * a) // (dr - d),
    (e * dr^2 - e * dr * d + rR * dr * a + dr * d * g - dr * r * a - d^2 * g) //
    (dr^2 + dr * a - d^2 - d * a),
]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

# QY system.
# (this is a big one)
ode = StructuralIdentifiability.@ODEmodel(
    P3'(t) = P4(t),
    P0'(t) = P1(t),
    P5'(t) =
        (
            -P0(t) * beta_SI * phi * Mar * Ks * siga2 + P0(t) * beta_SI * Mar * Ks * siga2 - P0(t) * phi * M * Mar * Ks * beta_SA +
            P0(t) * phi * M * Ks * siga2 * beta_SA +
            P0(t) * M * Mar * Ks * beta_SA - P1(t) * beta_SI * phi * Mar * siga2 -
            P1(t) * beta_SI * phi * Ks * siga2 +
            P1(t) * beta_SI * Mar * siga2 +
            P1(t) * beta_SI * Ks * siga2 - P1(t) * phi * M * Mar * beta_SA +
            P1(t) * phi * M * siga2 * beta_SA - P1(t) * phi * Mar * Ks * beta_SA +
            P1(t) * phi * Ks * siga2 * beta_SA +
            P1(t) * M * Mar * beta_SA +
            P1(t) * M * Ks * beta_SA +
            P1(t) * Mar * Ks * beta_SA - beta_SI * phi * P2(t) * siga2 +
            beta_SI * P2(t) * siga2 +
            P3(t) * beta_SA - phi * M * Mar * P5(t) * siga2 - phi * M * beta * siga2 -
            phi * P2(t) * Mar * beta_SA +
            phi * P2(t) * siga2 * beta_SA +
            M * P2(t) * beta_SA +
            M * Mar * P5(t) * siga2 +
            M * beta * siga2 +
            P2(t) * Mar * beta_SA +
            P2(t) * Ks * beta_SA
        ) // (phi * M * siga2 - M * siga2),
    P4'(t) =
        (
            -siga1 * P0(t)^2 * beta_SI * phi * M * Mar * Ks^2 * siga2^2 +
            siga1 * P0(t)^2 * beta_SI * M * Mar * Ks^2 * siga2^2 -
            siga1 * P0(t)^2 * phi * M^2 * Mar * Ks^2 * siga2 * beta_SA +
            siga1 * P0(t)^2 * phi * M^2 * Ks^2 * siga2^2 * beta_SA +
            siga1 * P0(t)^2 * M^2 * Mar * Ks^2 * siga2 * beta_SA -
            siga1 * P0(t) * P1(t) * beta_SI * phi * M * Mar * Ks^2 * siga2 -
            2 * siga1 * P0(t) * P1(t) * beta_SI * phi * M * Mar * Ks * siga2^2 -
            siga1 * P0(t) * P1(t) * beta_SI * phi * M * Ks^2 * siga2^2 -
            siga1 * P0(t) * P1(t) * beta_SI * phi * Mar * Ks^2 * siga2^2 +
            siga1 * P0(t) * P1(t) * beta_SI * M * Mar * Ks^2 * siga2 +
            2 * siga1 * P0(t) * P1(t) * beta_SI * M * Mar * Ks * siga2^2 +
            siga1 * P0(t) * P1(t) * beta_SI * M * Ks^2 * siga2^2 +
            siga1 * P0(t) * P1(t) * beta_SI * Mar * Ks^2 * siga2^2 -
            siga1 * P0(t) * P1(t) * phi * M^2 * Mar * Ks^2 * beta_SA -
            2 * siga1 * P0(t) * P1(t) * phi * M^2 * Mar * Ks * siga2 * beta_SA +
            siga1 * P0(t) * P1(t) * phi * M^2 * Ks^2 * siga2 * beta_SA +
            2 * siga1 * P0(t) * P1(t) * phi * M^2 * Ks * siga2^2 * beta_SA -
            2 * siga1 * P0(t) * P1(t) * phi * M * Mar * Ks^2 * siga2 * beta_SA +
            2 * siga1 * P0(t) * P1(t) * phi * M * Ks^2 * siga2^2 * beta_SA +
            siga1 * P0(t) * P1(t) * M^2 * Mar * Ks^2 * beta_SA +
            2 * siga1 * P0(t) * P1(t) * M^2 * Mar * Ks * siga2 * beta_SA +
            siga1 * P0(t) * P1(t) * M^2 * Ks^2 * siga2 * beta_SA +
            2 * siga1 * P0(t) * P1(t) * M * Mar * Ks^2 * siga2 * beta_SA -
            siga1 * P0(t) * beta_SI * P3(t) * phi * Mar * Ks * siga2 +
            siga1 * P0(t) * beta_SI * P3(t) * Mar * Ks * siga2 -
            siga1 * P0(t) * beta_SI * phi * M * P2(t) * Mar * Ks * siga2 -
            siga1 * P0(t) * beta_SI * phi * M * P2(t) * Ks * siga2^2 -
            siga1 * P0(t) * beta_SI * phi * P2(t) * Mar * Ks^2 * siga2 -
            siga1 * P0(t) * beta_SI * phi * P2(t) * Mar * Ks * siga2^2 +
            siga1 * P0(t) * beta_SI * M * P2(t) * Mar * Ks * siga2 +
            siga1 * P0(t) * beta_SI * M * P2(t) * Ks * siga2^2 +
            siga1 * P0(t) * beta_SI * P2(t) * Mar * Ks^2 * siga2 +
            siga1 * P0(t) * beta_SI * P2(t) * Mar * Ks * siga2^2 -
            siga1 * P0(t) * P3(t) * phi * M * Mar * Ks * beta_SA +
            siga1 * P0(t) * P3(t) * phi * M * Ks * siga2 * beta_SA +
            siga1 * P0(t) * P3(t) * M * Mar * Ks * beta_SA +
            siga1 * P0(t) * P3(t) * M * Ks * siga2 * beta_SA -
            siga1 * P0(t) * phi * M^2 * P2(t) * Mar * Ks * beta_SA +
            siga1 * P0(t) * phi * M^2 * P2(t) * Ks * siga2 * beta_SA -
            siga1 * P0(t) * phi * M^2 * Mar * P5(t) * Ks * siga2^2 -
            siga1 * P0(t) * phi * M^2 * Ks * beta * siga2^2 -
            siga1 * P0(t) * phi * M * P2(t) * Mar * Ks^2 * beta_SA -
            2 * siga1 * P0(t) * phi * M * P2(t) * Mar * Ks * siga2 * beta_SA +
            siga1 * P0(t) * phi * M * P2(t) * Ks^2 * siga2 * beta_SA +
            2 * siga1 * P0(t) * phi * M * P2(t) * Ks * siga2^2 * beta_SA +
            siga1 * P0(t) * M^2 * P2(t) * Mar * Ks * beta_SA +
            siga1 * P0(t) * M^2 * P2(t) * Ks * siga2 * beta_SA +
            siga1 * P0(t) * M^2 * Mar * P5(t) * Ks * siga2^2 +
            siga1 * P0(t) * M^2 * Ks * beta * siga2^2 +
            siga1 * P0(t) * M * P2(t) * Mar * Ks^2 * beta_SA +
            2 * siga1 * P0(t) * M * P2(t) * Mar * Ks * siga2 * beta_SA +
            siga1 * P0(t) * M * P2(t) * Ks^2 * siga2 * beta_SA -
            siga1 * P1(t)^2 * beta_SI * phi * M * Mar * Ks * siga2 -
            siga1 * P1(t)^2 * beta_SI * phi * M * Mar * siga2^2 -
            siga1 * P1(t)^2 * beta_SI * phi * M * Ks^2 * siga2 -
            siga1 * P1(t)^2 * beta_SI * phi * M * Ks * siga2^2 -
            siga1 * P1(t)^2 * beta_SI * phi * Mar * Ks * siga2^2 -
            siga1 * P1(t)^2 * beta_SI * phi * Ks^2 * siga2^2 +
            siga1 * P1(t)^2 * beta_SI * M * Mar * Ks * siga2 +
            siga1 * P1(t)^2 * beta_SI * M * Mar * siga2^2 +
            siga1 * P1(t)^2 * beta_SI * M * Ks^2 * siga2 +
            siga1 * P1(t)^2 * beta_SI * M * Ks * siga2^2 +
            siga1 * P1(t)^2 * beta_SI * Mar * Ks * siga2^2 +
            siga1 * P1(t)^2 * beta_SI * Ks^2 * siga2^2 -
            siga1 * P1(t)^2 * phi * M^2 * Mar * Ks * beta_SA -
            siga1 * P1(t)^2 * phi * M^2 * Mar * siga2 * beta_SA +
            siga1 * P1(t)^2 * phi * M^2 * Ks * siga2 * beta_SA +
            siga1 * P1(t)^2 * phi * M^2 * siga2^2 * beta_SA -
            siga1 * P1(t)^2 * phi * M * Mar * Ks^2 * beta_SA -
            2 * siga1 * P1(t)^2 * phi * M * Mar * Ks * siga2 * beta_SA +
            siga1 * P1(t)^2 * phi * M * Ks^2 * siga2 * beta_SA +
            2 * siga1 * P1(t)^2 * phi * M * Ks * siga2^2 * beta_SA -
            siga1 * P1(t)^2 * phi * Mar * Ks^2 * siga2 * beta_SA +
            siga1 * P1(t)^2 * phi * Ks^2 * siga2^2 * beta_SA +
            siga1 * P1(t)^2 * M^2 * Mar * Ks * beta_SA +
            siga1 * P1(t)^2 * M^2 * Mar * siga2 * beta_SA +
            siga1 * P1(t)^2 * M^2 * Ks^2 * beta_SA +
            siga1 * P1(t)^2 * M^2 * Ks * siga2 * beta_SA +
            siga1 * P1(t)^2 * M * Mar * Ks^2 * beta_SA +
            2 * siga1 * P1(t)^2 * M * Mar * Ks * siga2 * beta_SA +
            siga1 * P1(t)^2 * M * Ks^2 * siga2 * beta_SA +
            siga1 * P1(t)^2 * Mar * Ks^2 * siga2 * beta_SA -
            siga1 * P1(t) * beta_SI * P3(t) * phi * Mar * siga2 -
            siga1 * P1(t) * beta_SI * P3(t) * phi * Ks * siga2 +
            siga1 * P1(t) * beta_SI * P3(t) * Mar * siga2 +
            siga1 * P1(t) * beta_SI * P3(t) * Ks * siga2 -
            siga1 * P1(t) * beta_SI * phi * M * P2(t) * Mar * siga2 -
            2 * siga1 * P1(t) * beta_SI * phi * M * P2(t) * Ks * siga2 -
            siga1 * P1(t) * beta_SI * phi * M * P2(t) * siga2^2 -
            siga1 * P1(t) * beta_SI * phi * P2(t) * Mar * Ks * siga2 -
            siga1 * P1(t) * beta_SI * phi * P2(t) * Mar * siga2^2 -
            siga1 * P1(t) * beta_SI * phi * P2(t) * Ks^2 * siga2 -
            2 * siga1 * P1(t) * beta_SI * phi * P2(t) * Ks * siga2^2 +
            siga1 * P1(t) * beta_SI * M * P2(t) * Mar * siga2 +
            2 * siga1 * P1(t) * beta_SI * M * P2(t) * Ks * siga2 +
            siga1 * P1(t) * beta_SI * M * P2(t) * siga2^2 +
            siga1 * P1(t) * beta_SI * P2(t) * Mar * Ks * siga2 +
            siga1 * P1(t) * beta_SI * P2(t) * Mar * siga2^2 +
            siga1 * P1(t) * beta_SI * P2(t) * Ks^2 * siga2 +
            2 * siga1 * P1(t) * beta_SI * P2(t) * Ks * siga2^2 -
            siga1 * P1(t) * P3(t) * phi * M * Mar * beta_SA +
            siga1 * P1(t) * P3(t) * phi * M * siga2 * beta_SA -
            siga1 * P1(t) * P3(t) * phi * Mar * Ks * beta_SA +
            siga1 * P1(t) * P3(t) * phi * Ks * siga2 * beta_SA +
            siga1 * P1(t) * P3(t) * M * Mar * beta_SA +
            2 * siga1 * P1(t) * P3(t) * M * Ks * beta_SA +
            siga1 * P1(t) * P3(t) * M * siga2 * beta_SA +
            siga1 * P1(t) * P3(t) * Mar * Ks * beta_SA +
            siga1 * P1(t) * P3(t) * Ks * siga2 * beta_SA -
            siga1 * P1(t) * phi * M^2 * P2(t) * Mar * beta_SA +
            siga1 * P1(t) * phi * M^2 * P2(t) * siga2 * beta_SA -
            siga1 * P1(t) * phi * M^2 * Mar * P5(t) * Ks * siga2 -
            siga1 * P1(t) * phi * M^2 * Mar * P5(t) * siga2^2 -
            siga1 * P1(t) * phi * M^2 * Ks * beta * siga2 -
            siga1 * P1(t) * phi * M^2 * Ks * siga2^2 -
            siga1 * P1(t) * phi * M^2 * beta * siga2^2 -
            3 * siga1 * P1(t) * phi * M * P2(t) * Mar * Ks * beta_SA -
            2 * siga1 * P1(t) * phi * M * P2(t) * Mar * siga2 * beta_SA +
            3 * siga1 * P1(t) * phi * M * P2(t) * Ks * siga2 * beta_SA +
            2 * siga1 * P1(t) * phi * M * P2(t) * siga2^2 * beta_SA -
            siga1 * P1(t) * phi * M * Mar * P5(t) * Ks * siga2^2 -
            siga1 * P1(t) * phi * M * Ks * beta * siga2^2 -
            siga1 * P1(t) * phi * P2(t) * Mar * Ks^2 * beta_SA -
            2 * siga1 * P1(t) * phi * P2(t) * Mar * Ks * siga2 * beta_SA +
            siga1 * P1(t) * phi * P2(t) * Ks^2 * siga2 * beta_SA +
            2 * siga1 * P1(t) * phi * P2(t) * Ks * siga2^2 * beta_SA +
            siga1 * P1(t) * M^2 * P2(t) * Mar * beta_SA +
            2 * siga1 * P1(t) * M^2 * P2(t) * Ks * beta_SA +
            siga1 * P1(t) * M^2 * P2(t) * siga2 * beta_SA +
            siga1 * P1(t) * M^2 * Mar * P5(t) * Ks * siga2 +
            siga1 * P1(t) * M^2 * Mar * P5(t) * siga2^2 +
            siga1 * P1(t) * M^2 * Ks * beta * siga2 +
            siga1 * P1(t) * M^2 * Ks * siga2^2 +
            siga1 * P1(t) * M^2 * beta * siga2^2 +
            3 * siga1 * P1(t) * M * P2(t) * Mar * Ks * beta_SA +
            2 * siga1 * P1(t) * M * P2(t) * Mar * siga2 * beta_SA +
            2 * siga1 * P1(t) * M * P2(t) * Ks^2 * beta_SA +
            3 * siga1 * P1(t) * M * P2(t) * Ks * siga2 * beta_SA +
            siga1 * P1(t) * M * Mar * P5(t) * Ks * siga2^2 +
            siga1 * P1(t) * M * Ks * beta * siga2^2 +
            siga1 * P1(t) * P2(t) * Mar * Ks^2 * beta_SA +
            2 * siga1 * P1(t) * P2(t) * Mar * Ks * siga2 * beta_SA +
            siga1 * P1(t) * P2(t) * Ks^2 * siga2 * beta_SA -
            siga1 * beta_SI * P3(t) * phi * P2(t) * siga2 +
            siga1 * beta_SI * P3(t) * P2(t) * siga2 -
            siga1 * beta_SI * phi * M * P2(t)^2 * siga2 -
            siga1 * beta_SI * phi * P2(t)^2 * Ks * siga2 -
            siga1 * beta_SI * phi * P2(t)^2 * siga2^2 +
            siga1 * beta_SI * M * P2(t)^2 * siga2 +
            siga1 * beta_SI * P2(t)^2 * Ks * siga2 +
            siga1 * beta_SI * P2(t)^2 * siga2^2 +
            siga1 * P3(t)^2 * beta_SA - siga1 * P3(t) * phi * M^2 * siga2 -
            siga1 * P3(t) * phi * M * Mar * P5(t) * siga2 -
            siga1 * P3(t) * phi * M * Ks * siga2 -
            siga1 * P3(t) * phi * M * beta * siga2 - siga1 * P3(t) * phi * M * siga2^2 -
            siga1 * P3(t) * phi * P2(t) * Mar * beta_SA +
            siga1 * P3(t) * phi * P2(t) * siga2 * beta_SA +
            siga1 * P3(t) * M^2 * siga2 +
            2 * siga1 * P3(t) * M * P2(t) * beta_SA +
            siga1 * P3(t) * M * Mar * P5(t) * siga2 +
            siga1 * P3(t) * M * Ks * siga2 +
            siga1 * P3(t) * M * beta * siga2 +
            siga1 * P3(t) * M * siga2^2 +
            siga1 * P3(t) * P2(t) * Mar * beta_SA +
            2 * siga1 * P3(t) * P2(t) * Ks * beta_SA +
            siga1 * P3(t) * P2(t) * siga2 * beta_SA -
            siga1 * phi * M^2 * P2(t) * Mar * P5(t) * siga2 -
            siga1 * phi * M^2 * P2(t) * Ks * siga2 -
            siga1 * phi * M^2 * P2(t) * beta * siga2 -
            siga1 * phi * M^2 * P2(t) * siga2^2 - siga1 * phi * M * P4(t) * siga2 -
            siga1 * phi * M * P2(t)^2 * Mar * beta_SA +
            siga1 * phi * M * P2(t)^2 * siga2 * beta_SA -
            siga1 * phi * M * P2(t) * Mar * P5(t) * Ks * siga2 -
            siga1 * phi * M * P2(t) * Mar * P5(t) * siga2^2 -
            siga1 * phi * M * P2(t) * Ks * beta * siga2 -
            siga1 * phi * M * P2(t) * Ks * siga2^2 -
            siga1 * phi * M * P2(t) * beta * siga2^2 -
            siga1 * phi * P2(t)^2 * Mar * Ks * beta_SA -
            siga1 * phi * P2(t)^2 * Mar * siga2 * beta_SA +
            siga1 * phi * P2(t)^2 * Ks * siga2 * beta_SA +
            siga1 * phi * P2(t)^2 * siga2^2 * beta_SA +
            siga1 * M^2 * P2(t)^2 * beta_SA +
            siga1 * M^2 * P2(t) * Mar * P5(t) * siga2 +
            siga1 * M^2 * P2(t) * Ks * siga2 +
            siga1 * M^2 * P2(t) * beta * siga2 +
            siga1 * M^2 * P2(t) * siga2^2 +
            siga1 * M * P4(t) * siga2 +
            siga1 * M * P2(t)^2 * Mar * beta_SA +
            2 * siga1 * M * P2(t)^2 * Ks * beta_SA +
            siga1 * M * P2(t)^2 * siga2 * beta_SA +
            siga1 * M * P2(t) * Mar * P5(t) * Ks * siga2 +
            siga1 * M * P2(t) * Mar * P5(t) * siga2^2 +
            siga1 * M * P2(t) * Ks * beta * siga2 +
            siga1 * M * P2(t) * Ks * siga2^2 +
            siga1 * M * P2(t) * beta * siga2^2 +
            siga1 * P2(t)^2 * Mar * Ks * beta_SA +
            siga1 * P2(t)^2 * Mar * siga2 * beta_SA +
            siga1 * P2(t)^2 * Ks^2 * beta_SA +
            siga1 * P2(t)^2 * Ks * siga2 * beta_SA -
            P0(t) * P1(t) * beta_SI * phi * M * Mar * Ks^2 * siga2^2 +
            P0(t) * P1(t) * beta_SI * M * Mar * Ks^2 * siga2^2 -
            P0(t) * P1(t) * phi * M^2 * Mar * Ks^2 * siga2 * beta_SA +
            P0(t) * P1(t) * phi * M^2 * Ks^2 * siga2^2 * beta_SA +
            P0(t) * P1(t) * M^2 * Mar * Ks^2 * siga2 * beta_SA -
            P0(t) * beta_SI * P3(t) * phi * M * Mar * Ks * siga2 -
            P0(t) * beta_SI * P3(t) * phi * Mar * Ks^2 * siga2 -
            P0(t) * beta_SI * P3(t) * phi * Mar * Ks * siga2^2 +
            P0(t) * beta_SI * P3(t) * M * Mar * Ks * siga2 +
            P0(t) * beta_SI * P3(t) * Mar * Ks^2 * siga2 +
            P0(t) * beta_SI * P3(t) * Mar * Ks * siga2^2 -
            P0(t) * beta_SI * phi * alpa * Mar * Ks * siga2 -
            P0(t) * beta_SI * phi * M * P2(t) * Mar * Ks^2 * siga2 -
            P0(t) * beta_SI * phi * M * P2(t) * Mar * Ks * siga2^2 -
            P0(t) * beta_SI * phi * P4(t) * Mar * Ks * siga2 -
            P0(t) * beta_SI * phi * P2(t) * Mar * Ks^2 * siga2^2 +
            P0(t) * beta_SI * alpa * Mar * Ks * siga2 +
            P0(t) * beta_SI * M * P2(t) * Mar * Ks^2 * siga2 +
            P0(t) * beta_SI * M * P2(t) * Mar * Ks * siga2^2 +
            P0(t) * beta_SI * P4(t) * Mar * Ks * siga2 +
            P0(t) * beta_SI * P2(t) * Mar * Ks^2 * siga2^2 -
            P0(t) * P3(t) * phi * M^2 * Mar * Ks * beta_SA +
            P0(t) * P3(t) * phi * M^2 * Ks * siga2 * beta_SA -
            P0(t) * P3(t) * phi * M * Mar * Ks^2 * beta_SA -
            P0(t) * P3(t) * phi * M * Mar * Ks * siga2 * beta_SA +
            P0(t) * P3(t) * phi * M * Ks^2 * siga2 * beta_SA +
            P0(t) * P3(t) * phi * M * Ks * siga2^2 * beta_SA +
            P0(t) * P3(t) * M^2 * Mar * Ks * beta_SA +
            P0(t) * P3(t) * M * Mar * Ks^2 * beta_SA +
            P0(t) * P3(t) * M * Mar * Ks * siga2 * beta_SA -
            P0(t) * phi * alpa * M * Mar * Ks * beta_SA +
            P0(t) * phi * alpa * M * Ks * siga2 * beta_SA -
            P0(t) * phi * M^2 * P2(t) * Mar * Ks^2 * beta_SA -
            P0(t) * phi * M^2 * P2(t) * Mar * Ks * siga2 * beta_SA +
            P0(t) * phi * M^2 * P2(t) * Ks^2 * siga2 * beta_SA +
            P0(t) * phi * M^2 * P2(t) * Ks * siga2^2 * beta_SA -
            P0(t) * phi * M * P4(t) * Mar * Ks * beta_SA +
            P0(t) * phi * M * P4(t) * Ks * siga2 * beta_SA -
            P0(t) * phi * M * P2(t) * Mar * Ks^2 * siga2 * beta_SA +
            P0(t) * phi * M * P2(t) * Ks^2 * siga2^2 * beta_SA +
            P0(t) * alpa * M * Mar * Ks * beta_SA +
            P0(t) * M^2 * P2(t) * Mar * Ks^2 * beta_SA +
            P0(t) * M^2 * P2(t) * Mar * Ks * siga2 * beta_SA +
            P0(t) * M * P4(t) * Mar * Ks * beta_SA +
            P0(t) * M * P2(t) * Mar * Ks^2 * siga2 * beta_SA -
            P1(t)^2 * beta_SI * phi * M * Mar * Ks * siga2^2 -
            P1(t)^2 * beta_SI * phi * M * Ks^2 * siga2^2 +
            P1(t)^2 * beta_SI * M * Mar * Ks * siga2^2 +
            P1(t)^2 * beta_SI * M * Ks^2 * siga2^2 -
            P1(t)^2 * phi * M^2 * Mar * Ks * siga2 * beta_SA +
            P1(t)^2 * phi * M^2 * Ks * siga2^2 * beta_SA -
            P1(t)^2 * phi * M * Mar * Ks^2 * siga2 * beta_SA +
            P1(t)^2 * phi * M * Ks^2 * siga2^2 * beta_SA +
            P1(t)^2 * M^2 * Mar * Ks * siga2 * beta_SA +
            P1(t)^2 * M^2 * Ks^2 * siga2 * beta_SA +
            P1(t)^2 * M * Mar * Ks^2 * siga2 * beta_SA -
            P1(t) * beta_SI * P3(t) * phi * M * Mar * siga2 -
            P1(t) * beta_SI * P3(t) * phi * M * Ks * siga2 -
            P1(t) * beta_SI * P3(t) * phi * Mar * Ks * siga2 -
            P1(t) * beta_SI * P3(t) * phi * Mar * siga2^2 -
            P1(t) * beta_SI * P3(t) * phi * Ks^2 * siga2 -
            P1(t) * beta_SI * P3(t) * phi * Ks * siga2^2 +
            P1(t) * beta_SI * P3(t) * M * Mar * siga2 +
            P1(t) * beta_SI * P3(t) * M * Ks * siga2 +
            P1(t) * beta_SI * P3(t) * Mar * Ks * siga2 +
            P1(t) * beta_SI * P3(t) * Mar * siga2^2 +
            P1(t) * beta_SI * P3(t) * Ks^2 * siga2 +
            P1(t) * beta_SI * P3(t) * Ks * siga2^2 -
            P1(t) * beta_SI * phi * alpa * Mar * siga2 -
            P1(t) * beta_SI * phi * alpa * Ks * siga2 -
            P1(t) * beta_SI * phi * M * P2(t) * Mar * Ks * siga2 -
            P1(t) * beta_SI * phi * M * P2(t) * Mar * siga2^2 -
            P1(t) * beta_SI * phi * M * P2(t) * Ks^2 * siga2 -
            2 * P1(t) * beta_SI * phi * M * P2(t) * Ks * siga2^2 -
            P1(t) * beta_SI * phi * P4(t) * Mar * siga2 -
            P1(t) * beta_SI * phi * P4(t) * Ks * siga2 -
            P1(t) * beta_SI * phi * P2(t) * Mar * Ks * siga2^2 -
            P1(t) * beta_SI * phi * P2(t) * Ks^2 * siga2^2 +
            P1(t) * beta_SI * alpa * Mar * siga2 +
            P1(t) * beta_SI * alpa * Ks * siga2 +
            P1(t) * beta_SI * M * P2(t) * Mar * Ks * siga2 +
            P1(t) * beta_SI * M * P2(t) * Mar * siga2^2 +
            P1(t) * beta_SI * M * P2(t) * Ks^2 * siga2 +
            2 * P1(t) * beta_SI * M * P2(t) * Ks * siga2^2 +
            P1(t) * beta_SI * P4(t) * Mar * siga2 +
            P1(t) * beta_SI * P4(t) * Ks * siga2 +
            P1(t) * beta_SI * P2(t) * Mar * Ks * siga2^2 +
            P1(t) * beta_SI * P2(t) * Ks^2 * siga2^2 -
            P1(t) * P3(t) * phi * M^2 * Mar * beta_SA +
            P1(t) * P3(t) * phi * M^2 * siga2 * beta_SA -
            2 * P1(t) * P3(t) * phi * M * Mar * Ks * beta_SA -
            P1(t) * P3(t) * phi * M * Mar * siga2 * beta_SA +
            2 * P1(t) * P3(t) * phi * M * Ks * siga2 * beta_SA +
            P1(t) * P3(t) * phi * M * siga2^2 * beta_SA -
            P1(t) * P3(t) * phi * Mar * Ks^2 * beta_SA -
            P1(t) * P3(t) * phi * Mar * Ks * siga2 * beta_SA +
            P1(t) * P3(t) * phi * Ks^2 * siga2 * beta_SA +
            P1(t) * P3(t) * phi * Ks * siga2^2 * beta_SA +
            P1(t) * P3(t) * M^2 * Mar * beta_SA +
            P1(t) * P3(t) * M^2 * Ks * beta_SA +
            2 * P1(t) * P3(t) * M * Mar * Ks * beta_SA +
            P1(t) * P3(t) * M * Mar * siga2 * beta_SA +
            P1(t) * P3(t) * M * Ks^2 * beta_SA +
            2 * P1(t) * P3(t) * M * Ks * siga2 * beta_SA +
            P1(t) * P3(t) * Mar * Ks^2 * beta_SA +
            P1(t) * P3(t) * Mar * Ks * siga2 * beta_SA -
            P1(t) * phi * alpa * M * Mar * beta_SA +
            P1(t) * phi * alpa * M * siga2 * beta_SA -
            P1(t) * phi * alpa * Mar * Ks * beta_SA +
            P1(t) * phi * alpa * Ks * siga2 * beta_SA -
            P1(t) * phi * M^2 * P2(t) * Mar * Ks * beta_SA -
            P1(t) * phi * M^2 * P2(t) * Mar * siga2 * beta_SA +
            P1(t) * phi * M^2 * P2(t) * Ks * siga2 * beta_SA +
            P1(t) * phi * M^2 * P2(t) * siga2^2 * beta_SA -
            P1(t) * phi * M^2 * Mar * P5(t) * Ks * siga2^2 -
            P1(t) * phi * M^2 * Ks * beta * siga2^2 -
            P1(t) * phi * M * P4(t) * Mar * beta_SA +
            P1(t) * phi * M * P4(t) * siga2 * beta_SA -
            P1(t) * phi * M * P2(t) * Mar * Ks^2 * beta_SA -
            3 * P1(t) * phi * M * P2(t) * Mar * Ks * siga2 * beta_SA +
            P1(t) * phi * M * P2(t) * Ks^2 * siga2 * beta_SA +
            3 * P1(t) * phi * M * P2(t) * Ks * siga2^2 * beta_SA -
            P1(t) * phi * P4(t) * Mar * Ks * beta_SA +
            P1(t) * phi * P4(t) * Ks * siga2 * beta_SA -
            P1(t) * phi * P2(t) * Mar * Ks^2 * siga2 * beta_SA +
            P1(t) * phi * P2(t) * Ks^2 * siga2^2 * beta_SA +
            P1(t) * alpa * M * Mar * beta_SA +
            P1(t) * alpa * M * Ks * beta_SA +
            P1(t) * alpa * Mar * Ks * beta_SA +
            P1(t) * M^2 * P2(t) * Mar * Ks * beta_SA +
            P1(t) * M^2 * P2(t) * Mar * siga2 * beta_SA +
            P1(t) * M^2 * P2(t) * Ks^2 * beta_SA +
            2 * P1(t) * M^2 * P2(t) * Ks * siga2 * beta_SA +
            P1(t) * M^2 * Mar * P5(t) * Ks * siga2^2 +
            P1(t) * M^2 * Ks * beta * siga2^2 +
            P1(t) * M * P4(t) * Mar * beta_SA +
            P1(t) * M * P4(t) * Ks * beta_SA +
            P1(t) * M * P2(t) * Mar * Ks^2 * beta_SA +
            3 * P1(t) * M * P2(t) * Mar * Ks * siga2 * beta_SA +
            2 * P1(t) * M * P2(t) * Ks^2 * siga2 * beta_SA +
            P1(t) * P4(t) * Mar * Ks * beta_SA +
            P1(t) * P2(t) * Mar * Ks^2 * siga2 * beta_SA -
            beta_SI * P3(t) * phi * M * P2(t) * siga2 -
            beta_SI * P3(t) * phi * P2(t) * Ks * siga2 -
            beta_SI * P3(t) * phi * P2(t) * siga2^2 +
            beta_SI * P3(t) * M * P2(t) * siga2 +
            beta_SI * P3(t) * P2(t) * Ks * siga2 +
            beta_SI * P3(t) * P2(t) * siga2^2 - beta_SI * phi * alpa * P2(t) * siga2 -
            beta_SI * phi * M * P2(t)^2 * Ks * siga2 -
            beta_SI * phi * M * P2(t)^2 * siga2^2 -
            beta_SI * phi * P4(t) * P2(t) * siga2 -
            beta_SI * phi * P2(t)^2 * Ks * siga2^2 +
            beta_SI * alpa * P2(t) * siga2 +
            beta_SI * M * P2(t)^2 * Ks * siga2 +
            beta_SI * M * P2(t)^2 * siga2^2 +
            beta_SI * P4(t) * P2(t) * siga2 +
            beta_SI * P2(t)^2 * Ks * siga2^2 +
            P3(t)^2 * M * beta_SA +
            P3(t)^2 * Ks * beta_SA +
            P3(t)^2 * siga2 * beta_SA - P3(t) * phi * M^2 * Mar * P5(t) * siga2 -
            P3(t) * phi * M^2 * Ks * siga2 - P3(t) * phi * M^2 * beta * siga2 -
            P3(t) * phi * M^2 * siga2^2 - P3(t) * phi * M * P2(t) * Mar * beta_SA +
            P3(t) * phi * M * P2(t) * siga2 * beta_SA -
            P3(t) * phi * M * Mar * P5(t) * Ks * siga2 -
            P3(t) * phi * M * Mar * P5(t) * siga2^2 -
            P3(t) * phi * M * Ks * beta * siga2 - P3(t) * phi * M * Ks * siga2^2 -
            P3(t) * phi * M * beta * siga2^2 -
            P3(t) * phi * P2(t) * Mar * Ks * beta_SA -
            P3(t) * phi * P2(t) * Mar * siga2 * beta_SA +
            P3(t) * phi * P2(t) * Ks * siga2 * beta_SA +
            P3(t) * phi * P2(t) * siga2^2 * beta_SA +
            P3(t) * alpa * beta_SA +
            P3(t) * M^2 * P2(t) * beta_SA +
            P3(t) * M^2 * Mar * P5(t) * siga2 +
            P3(t) * M^2 * Ks * siga2 +
            P3(t) * M^2 * beta * siga2 +
            P3(t) * M^2 * siga2^2 +
            P3(t) * M * P2(t) * Mar * beta_SA +
            3 * P3(t) * M * P2(t) * Ks * beta_SA +
            2 * P3(t) * M * P2(t) * siga2 * beta_SA +
            P3(t) * M * Mar * P5(t) * Ks * siga2 +
            P3(t) * M * Mar * P5(t) * siga2^2 +
            P3(t) * M * Ks * beta * siga2 +
            P3(t) * M * Ks * siga2^2 +
            P3(t) * M * beta * siga2^2 +
            P3(t) * P4(t) * beta_SA +
            P3(t) * P2(t) * Mar * Ks * beta_SA +
            P3(t) * P2(t) * Mar * siga2 * beta_SA +
            P3(t) * P2(t) * Ks^2 * beta_SA +
            2 * P3(t) * P2(t) * Ks * siga2 * beta_SA -
            phi * alpa * M * Mar * P5(t) * siga2 - phi * alpa * M * beta * siga2 -
            phi * alpa * P2(t) * Mar * beta_SA + phi * alpa * P2(t) * siga2 * beta_SA -
            phi * M^2 * P4(t) * siga2 - phi * M^2 * P2(t) * Mar * P5(t) * Ks * siga2 -
            phi * M^2 * P2(t) * Mar * P5(t) * siga2^2 -
            phi * M^2 * P2(t) * Ks * beta * siga2 - phi * M^2 * P2(t) * Ks * siga2^2 -
            phi * M^2 * P2(t) * beta * siga2^2 - phi * M * P4(t) * Mar * P5(t) * siga2 -
            phi * M * P4(t) * Ks * siga2 - phi * M * P4(t) * beta * siga2 -
            phi * M * P4(t) * siga2^2 - phi * M * P2(t)^2 * Mar * Ks * beta_SA -
            phi * M * P2(t)^2 * Mar * siga2 * beta_SA +
            phi * M * P2(t)^2 * Ks * siga2 * beta_SA +
            phi * M * P2(t)^2 * siga2^2 * beta_SA -
            phi * M * P2(t) * Mar * P5(t) * Ks * siga2^2 -
            phi * M * P2(t) * Ks * beta * siga2^2 -
            phi * P4(t) * P2(t) * Mar * beta_SA +
            phi * P4(t) * P2(t) * siga2 * beta_SA -
            phi * P2(t)^2 * Mar * Ks * siga2 * beta_SA +
            phi * P2(t)^2 * Ks * siga2^2 * beta_SA +
            alpa * M * P2(t) * beta_SA +
            alpa * M * Mar * P5(t) * siga2 +
            alpa * M * beta * siga2 +
            alpa * P2(t) * Mar * beta_SA +
            alpa * P2(t) * Ks * beta_SA +
            M^2 * P4(t) * siga2 +
            M^2 * P2(t)^2 * Ks * beta_SA +
            M^2 * P2(t)^2 * siga2 * beta_SA +
            M^2 * P2(t) * Mar * P5(t) * Ks * siga2 +
            M^2 * P2(t) * Mar * P5(t) * siga2^2 +
            M^2 * P2(t) * Ks * beta * siga2 +
            M^2 * P2(t) * Ks * siga2^2 +
            M^2 * P2(t) * beta * siga2^2 +
            M * P4(t) * P2(t) * beta_SA +
            M * P4(t) * Mar * P5(t) * siga2 +
            M * P4(t) * Ks * siga2 +
            M * P4(t) * beta * siga2 +
            M * P4(t) * siga2^2 +
            M * P2(t)^2 * Mar * Ks * beta_SA +
            M * P2(t)^2 * Mar * siga2 * beta_SA +
            M * P2(t)^2 * Ks^2 * beta_SA +
            2 * M * P2(t)^2 * Ks * siga2 * beta_SA +
            M * P2(t) * Mar * P5(t) * Ks * siga2^2 +
            M * P2(t) * Ks * beta * siga2^2 +
            P4(t) * P2(t) * Mar * beta_SA +
            P4(t) * P2(t) * Ks * beta_SA +
            P2(t)^2 * Mar * Ks * siga2 * beta_SA +
            P2(t)^2 * Ks^2 * siga2 * beta_SA
        ) // (phi * M * siga2 - M * siga2),
    P1'(t) = P2(t),
    P2'(t) = P3(t),
    y(t) = P0(t)
)
ident_funcs = [
    Mar,
    Ks,
    alpa,
    (siga1 + phi * Mar - Mar) // phi,
    (siga1 * phi - siga1 - phi * Mar + Mar) // (siga1 * phi * beta_SA),
    (
        siga1 * beta_SA + beta_SI * phi * siga2 + phi * Mar * beta_SA -
        phi * siga2 * beta_SA - Mar * beta_SA
    ) // (phi * M * siga2),
    (
        siga1 * phi * M * siga2 - siga1 * M * siga2 - phi * M * Mar * siga2 +
        M * Mar * siga2
    ) // (siga1 * beta_SA + phi * Mar * beta_SA - Mar * beta_SA),
    (
        siga1 * beta_SI * phi * siga2 - siga1 * beta_SI * siga2 -
        siga1 * phi * siga2 * beta_SA - siga1 * M * beta_SA -
        beta_SI * phi * Mar * siga2 + beta_SI * Mar * siga2 - phi * M * Mar * beta_SA +
        M * Mar * beta_SA
    ) // (siga1 * beta_SA + phi * Mar * beta_SA - Mar * beta_SA),
    (
        siga1^2 * beta_SA + siga1 * beta_SI * phi * siga2 - siga1 * beta_SI * siga2 +
        siga1 * phi * Mar * beta_SA - siga1 * phi * siga2 * beta_SA -
        siga1 * Mar * beta_SA +
        siga1 * siga2 * beta_SA +
        beta_SI * phi * M * siga2 - beta_SI * phi * Mar * siga2 +
        beta_SI * phi * siga2^2 +
        beta_SI * Mar * siga2 +
        phi * Mar * siga2 * beta_SA - phi * siga2^2 * beta_SA - Mar * siga2 * beta_SA
    ) // (phi * M * siga2),
]
# Really large and takes a lot of time, so commented
# push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

###
# Cases with states

ode = StructuralIdentifiability.@ODEmodel(x'(t) = x(t), y(t) = x(t))
T = typeof(x)
ident_funcs = [x]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs, with_states = true))

ode = StructuralIdentifiability.@ODEmodel(x'(t) = a * x(t) + u(t), y(t) = x(t))
ident_funcs = [a, x]
push!(test_cases, (ode = ode, with_states = true, ident_funcs = ident_funcs))

ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = a * x1(t) - b * x1(t) * x2(t),
    x2'(t) = -c * x2(t) + d * x1(t) * x2(t),
    y(t) = x1(t)
)
ident_funcs = [x1, c, d, a, b * x2]
push!(test_cases, (ode = ode, with_states = true, ident_funcs = ident_funcs))

ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = a * x1(t) - b * x1(t) * x2(t),
    x2'(t) = -c * x2(t) + d * x1(t) * x2(t),
    y(t) = x1(t)
)
ident_funcs = [x1, c, d, a, b * x2]
push!(test_cases, (ode = ode, with_states = true, ident_funcs = ident_funcs))

# Diagonal with simple spectrum and observable states
ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = λ1 * x1(t) + β1 * u1(t),
    x2'(t) = λ2 * x2(t) + β2 * u2(t),
    x3'(t) = λ3 * x3(t) + β3 * u3(t),
    y(t) = x1(t) + x2(t) + x3(t)
)
ident_funcs = [λ1, λ2, λ3, β1, β2, β3, x1, x2, x3]
push!(test_cases, (ode = ode, with_states = true, ident_funcs = ident_funcs))

ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = x1(t) + Θ * x2(t),
    x2'(t) = 0,
    y(t) = x1(t)
)
ident_funcs = [x1, x2 * Θ]
push!(test_cases, (ode = ode, with_states = true, ident_funcs = ident_funcs))

ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = α * x2(t),
    x2'(t) = x3(t),
    x3'(t) = C,
    y(t) = x1(t)
)
ident_funcs = [x1, α * x2, α * x3, α * C]
push!(test_cases, (ode = ode, with_states = true, ident_funcs = ident_funcs))

ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = α * (x1 - x2),
    x2'(t) = α * (x1 + x2),
    y(t) = (x1^2 + x2^2) // 2,
)
ident_funcs = [α, x1^2 + x2^2]
push!(test_cases, (ode = ode, with_states = true, ident_funcs = ident_funcs))

ode = StructuralIdentifiability.@ODEmodel(x'(t) = a * x(t) + b * u(t), y(t) = c * x(t))
ident_funcs = [a, b * c, x * c]
push!(test_cases, (ode = ode, with_states = true, ident_funcs = ident_funcs))

# llw1987 model
ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = -p1 * x1(t) + p2 * u(t),
    x2'(t) = -p3 * x2(t) + p4 * u(t),
    x3'(t) = -(p1 + p3) * x3(t) + (p4 * x1(t) + p2 * x2(t)) * u(t),
    y1(t) = x3(t)
)
ident_funcs = [
    x3,
    x2 * x1 // one(x1),
    p3 * p1 // one(x1),
    p2 * p4 // one(x1),
    (p3 + p1) // one(x1),
    (p2 * x2 + p4 * x1) // one(x1),
    (p2 * x2 - p4 * x1) // (p3 - p1),
]
push!(test_cases, (ode = ode, with_states = true, ident_funcs = ident_funcs))

# Regression test: Previously failed for with_states=true because of the bug in
# `linear_relations_between_normal_forms`
# Fujita
ode = StructuralIdentifiability.@ODEmodel(
    EGFR'(t) =
        EGFR_turnover * pro_EGFR(t) + EGF_EGFR(t) * reaction_1_k2 -
        EGFR(t) * EGFR_turnover - EGF_EGFR(t) * reaction_1_k1,
    pEGFR'(t) =
        EGF_EGFR(t) * reaction_9_k1 - pEGFR(t) * reaction_4_k1 +
        pEGFR_Akt(t) * reaction_2_k2 +
        pEGFR_Akt(t) * reaction_3_k1 - Akt(t) * pEGFR(t) * reaction_2_k1,
    pEGFR_Akt'(t) =
        Akt(t) * pEGFR(t) * reaction_2_k1 - pEGFR_Akt(t) * reaction_3_k1 -
        pEGFR_Akt(t) * reaction_2_k2,
    Akt'(t) =
        pAkt(t) * reaction_7_k1 + pEGFR_Akt(t) * reaction_2_k2 -
        Akt(t) * pEGFR(t) * reaction_2_k1,
    pAkt'(t) =
        pAkt_S6(t) * reaction_5_k2 - pAkt(t) * reaction_7_k1 +
        pAkt_S6(t) * reaction_6_k1 +
        pEGFR_Akt(t) * reaction_3_k1 - S6(t) * pAkt(t) * reaction_5_k1,
    S6'(t) =
        pAkt_S6(t) * reaction_5_k2 + pS6(t) * reaction_8_k1 -
        S6(t) * pAkt(t) * reaction_5_k1,
    pAkt_S6'(t) =
        S6(t) * pAkt(t) * reaction_5_k1 - pAkt_S6(t) * reaction_6_k1 -
        pAkt_S6(t) * reaction_5_k2,
    pS6'(t) = pAkt_S6(t) * reaction_6_k1 - pS6(t) * reaction_8_k1,
    EGF_EGFR'(t) =
        EGF_EGFR(t) * reaction_1_k1 - EGF_EGFR(t) * reaction_9_k1 -
        EGF_EGFR(t) * reaction_1_k2,
    y1(t) = a1 * (pEGFR(t) + pEGFR_Akt(t)),
    y2(t) = a2 * (pAkt(t) + pAkt_S6(t)),
    y3(t) = a3 * pS6(t)
)
ident_funcs = [
    (EGF_EGFR * reaction_9_k1) // pS6,
    reaction_8_k1,
    a3 // reaction_5_k1,
    reaction_3_k1,
    reaction_2_k2,
    reaction_2_k1 // reaction_5_k1,
    a1 // reaction_5_k1,
    pAkt * reaction_5_k1,
    pEGFR * reaction_5_k1,
    pS6 * reaction_5_k1,
    reaction_5_k2,
    reaction_6_k1,
    a2 // reaction_5_k1,
    reaction_7_k1,
    reaction_4_k1,
    pAkt_S6 * reaction_5_k1,
    reaction_1_k1 - reaction_1_k2 - reaction_9_k1,
    pEGFR_Akt * reaction_5_k1,
    S6 * reaction_5_k1,
    Akt * reaction_5_k1,
]
push!(test_cases, (ode = ode, with_states = true, ident_funcs = ident_funcs))

# Bruno2016 model 
ode = StructuralIdentifiability.@ODEmodel(
    beta'(t) = -kbeta * beta(t),
    cry'(t) = -cry(t) * kcrybeta - cry(t) * kcryOH,
    zea'(t) = -zea(t) * kzea,
    beta10'(t) = cry(t) * kcryOH - beta10(t) * kbeta10 + kbeta * beta(t),
    OHbeta10'(t) = cry(t) * kcrybeta + zea(t) * kzea - OHbeta10(t) * kOHbeta10,
    betaio'(t) = cry(t) * kcrybeta + beta10(t) * kbeta10 + kbeta * beta(t),
    OHbetaio'(t) = cry(t) * kcryOH + zea(t) * kzea + OHbeta10(t) * kOHbeta10,
    y1(t) = beta(t),
    y2(t) = beta10(t)
)
ident_funcs = [beta10, beta, kbeta, kbeta10, cry * kcryOH, kcrybeta + kcryOH]
push!(test_cases, (ode = ode, with_states = true, ident_funcs = ident_funcs))

# STAT-5 model from 
# MODELING THE NONLINEAR DYNAMICS OF CELLULAR SIGNAL TRANSDUCTION
# DOI: https://doi.org/10.1142/S0218127404010461
ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = -k1 * x1 * EpoR_A,
    x2'(t) = k1 * x1 * EpoR_A - k2 * x2^2,
    x3'(t) = -k3 * x3 + 0.5 * k2 * x2^2,
    x4'(t) = k3 * x3,
    y1(t) = k5 * (x2 + 2x3),
    y2(t) = k6 * (x1 + x2 + 2x3),
    y3(t) = k7 * EpoR_A
)
ident_funcs = [k2 // k6, k3, EpoR_A * k7, EpoR_A * k1, k5 // k6]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs))

ode = @ODEmodel(x1'(t) = x1, x2'(t) = x2, y(t) = x1 + x2(t))
ident_funcs = [x1 + x2]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs, with_states = true))

# SEUIR model from the reparametrization docs
# https://docs.sciml.ai/StructuralIdentifiability/stable/tutorials/reparametrization/
ode = @ODEmodel(
    S'(t) = -b * (U(t) + I(t)) * S(t) / N,
    E'(t) = b * (U(t) + I(t)) * S(t) / N - g * E(t),
    U'(t) = (1 - a) * g * E(t) - d * U(t),
    I'(t) = a * g * E(t) - d * I(t),
    y(t) = I(t)
)
ident_funcs = [I, d, g, S * a, E * a, U * a + I * a, (N * a) // b]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs, with_states = true))

ode = @ODEmodel(x'(t) = alpha * x(t)^2, y(t) = x(t)^2)
ident_funcs = [alpha^2, x * alpha]
push!(test_cases, (ode = ode, ident_funcs = ident_funcs, with_states = true))

# TODO: verify that Maple returns the same
@testset "Identifiable functions of parameters" begin
    p = 0.99
    for case in test_cases
        for simplify in [:weak, :standard] #:strong?
            ode = case.ode
            true_ident_funcs = case.ident_funcs
            with_states = false
            if haskey(case, :with_states)
                with_states = case.with_states
            end
            result_funcs = StructuralIdentifiability.find_identifiable_functions(
                ode,
                simplify = simplify,
                with_states = with_states,
            )

            if isempty(true_ident_funcs)
                @test isempty(result_funcs)
                continue
            end

            @test parent(numerator(result_funcs[1])) == parent(ode)

            R = parent(numerator(result_funcs[1]))

            @info "Test, result_funcs = \n$result_funcs" case simplify R with_states

            true_ident_funcs = map(f -> f // one(f), true_ident_funcs)
            true_ident_funcs = map(
                f -> StructuralIdentifiability.parent_ring_change(f, R),
                true_ident_funcs,
            )

            # Check inclusion <true funcs> in <result funcs>
            @test StructuralIdentifiability.fields_equal(
                StructuralIdentifiability.RationalFunctionField(result_funcs),
                StructuralIdentifiability.RationalFunctionField(true_ident_funcs),
                p,
            )
            if simplify != :weak
                @info gens(parent(numerator(first(result_funcs))))
                # To keep track of changes in the simplification
                @test Set(result_funcs) == Set(true_ident_funcs)
            end
        end
    end
end
