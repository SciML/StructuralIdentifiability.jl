using StructuralIdentifiability
using BenchmarkTools
import Nemo, Profile

macro myprof(ex)
    :((VSCodeServer.Profile).clear();
    Profile.init(n = 10^8, delay = 0.00001);
    Profile.@profile $ex;
    VSCodeServer.view_profile(;))
end

mapk_6_out = @ODEmodel(
    KS00'(t) =
        -a00 * K(t) * S00(t) +
        b00 * KS00(t) +
        gamma0100 * FS01(t) +
        gamma1000 * FS10(t) +
        gamma1100 * FS11(t),
    KS01'(t) =
        -a01 * K(t) * S01(t) + b01 * KS01(t) + c0001 * KS00(t) - alpha01 * F(t) * S01(t) +
        beta01 * FS01(t) +
        gamma1101 * FS11(t),
    KS10'(t) =
        -a10 * K(t) * S10(t) + b10 * KS10(t) + c0010 * KS00(t) - alpha10 * F(t) * S10(t) +
        beta10 * FS10(t) +
        gamma1110 * FS11(t),
    FS01'(t) =
        -alpha11 * F(t) * S11(t) +
        beta11 * FS11(t) +
        c0111 * KS01(t) +
        c1011 * KS10(t) +
        c0011 * KS00(t),
    FS10'(t) = a00 * K(t) * S00(t) - (b00 + c0001 + c0010 + c0011) * KS00(t),
    FS11'(t) = a01 * K(t) * S01(t) - (b01 + c0111) * KS01(t),
    K'(t) = a10 * K(t) * S10(t) - (b10 + c1011) * KS10(t),
    F'(t) = alpha01 * F(t) * S01(t) - (beta01 + gamma0100) * FS01(t),
    S00'(t) = alpha10 * F(t) * S10(t) - (beta10 + gamma1000) * FS10(t),
    S01'(t) =
        alpha11 * F(t) * S11(t) - (beta11 + gamma1101 + gamma1110 + gamma1100) * FS11(t),
    S10'(t) =
        -a00 * K(t) * S00(t) + (b00 + c0001 + c0010 + c0011) * KS00(t) -
        a01 * K(t) * S01(t) + (b01 + c0111) * KS01(t) - a10 * K(t) * S10(t) +
        (b10 + c1011) * KS10(t),
    S11'(t) =
        -alpha01 * F(t) * S01(t) + (beta01 + gamma0100) * FS01(t) -
        alpha10 * F(t) * S10(t) + (beta10 + gamma1000) * FS10(t) -
        alpha11 * F(t) * S11(t) +
        (beta11 + gamma1101 + gamma1110 + gamma1100) * FS11(t),
    y0(t) = K(t),
    y1(t) = F(t),
    y2(t) = S00(t),
    y3(t) = S01(t),
    y4(t) = S10(t),
    y5(t) = S11(t)
)

siwr = @ODEmodel(
    S'(t) = mu - bi * S(t) * I(t) - bw * S(t) * W(t) - mu * S(t) + a * R(t),
    I'(t) = bw * S(t) * W(t) + bi * S(t) * I(t) - (gam + mu) * I(t),
    W'(t) = xi * (I(t) - W(t)),
    R'(t) = gam * I(t) - (mu + a) * R(t),
    y(t) = k * I(t)
)

sirs = @ODEmodel(
    s'(t) = mu - mu * s(t) - b0 * (1 + b1 * x1(t)) * i(t) * s(t) + g * r(t),
    i'(t) = b0 * (1 + b1 * x1(t)) * i(t) * s(t) - (nu + mu) * i(t),
    r'(t) = nu * i(t) - (mu + g) * r(t),
    x1'(t) = -M * x2(t),
    x2'(t) = M * x1(t),
    y1(t) = i(t),
    y2(t) = r(t)
)

covid = @ODEmodel(
    S'(t) = -b * S(t) * (I(t) + J(t) + q * A(t)) * Ninv(t),
    E'(t) = b * S(t) * (I(t) + J(t) + q * A(t)) * Ninv(t) - k * E(t),
    A'(t) = k * (1 - r) * E(t) - g1 * A(t),
    I'(t) = k * r * E(t) - (alpha + g1) * I(t),
    J'(t) = alpha * I(t) - g2 * J(t),
    C'(t) = alpha * I(t),
    Ninv'(t) = 0,
    y(t) = C(t),
    y2(t) = Ninv(t)
)

pharm = @ODEmodel(
    x0'(t) =
        a1 * (x1(t) - x0(t)) - (ka * n * x0(t)) / (kc * ka + kc * x2(t) + ka * x0(t)),
    x1'(t) = a2 * (x0(t) - x1(t)),
    x2'(t) =
        b1 * (x3(t) - x2(t)) - (kc * n * x2(t)) / (kc * ka + kc * x2(t) + ka * x0(t)),
    x3'(t) = b2 * (x2(t) - x3(t)),
    y1(t) = x0(t)
)

St = @ODEmodel(
    S'(t) = r * S(t) - (e + a * W(t)) * S(t) - d * W(t) * S(t) + g * R(t),
    R'(t) = rR * R(t) + (e + a * W(t)) * S(t) - dr * W(t) * R(t) - g * R(t),
    W'(t) = Dd * (T - W(t)),
    y1(t) = S(t) + R(t),
    y2(t) = T
)

cd8 = @ODEmodel(
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

qwwc = @ODEmodel(
    x'(t) = a * (y(t) - x(t)) + y(t) * z(t),
    y'(t) = b * (x(t) + y(t)) - x(t) * z(t),
    z'(t) = -c * z(t) - d * w(t) + x(t) * y(t),
    w'(t) = e * z(t) - f * w(t) + x(t) * y(t),
    g(t) = x(t)
)

@time StructuralIdentifiability.find_identifiable_functions(qwwc)

@myprof find_identifiable_functions(qwwc)
