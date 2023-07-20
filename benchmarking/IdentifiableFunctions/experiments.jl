using StructuralIdentifiability
using BenchmarkTools
import Nemo, Profile

macro myprof(ex)
    :((VSCodeServer.Profile).clear();
    Profile.init(n = 10^8, delay = 0.00001);
    Profile.@profile $ex;
    VSCodeServer.view_profile(;))
end

mapk_6_out = StructuralIdentifiability.@ODEmodel(
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

siwr = StructuralIdentifiability.@ODEmodel(
    S'(t) = mu - bi * S(t) * I(t) - bw * S(t) * W(t) - mu * S(t) + a * R(t),
    I'(t) = bw * S(t) * W(t) + bi * S(t) * I(t) - (gam + mu) * I(t),
    W'(t) = xi * (I(t) - W(t)),
    R'(t) = gam * I(t) - (mu + a) * R(t),
    y(t) = k * I(t)
)

sirs = StructuralIdentifiability.@ODEmodel(
    s'(t) = mu - mu * s(t) - b0 * (1 + b1 * x1(t)) * i(t) * s(t) + g * r(t),
    i'(t) = b0 * (1 + b1 * x1(t)) * i(t) * s(t) - (nu + mu) * i(t),
    r'(t) = nu * i(t) - (mu + g) * r(t),
    x1'(t) = -M * x2(t),
    x2'(t) = M * x1(t),
    y1(t) = i(t),
    y2(t) = r(t)
)

covid = StructuralIdentifiability.@ODEmodel(
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

pharm = StructuralIdentifiability.@ODEmodel(
    x0'(t) =
        a1 * (x1(t) - x0(t)) - (ka * n * x0(t)) / (kc * ka + kc * x2(t) + ka * x0(t)),
    x1'(t) = a2 * (x0(t) - x1(t)),
    x2'(t) =
        b1 * (x3(t) - x2(t)) - (kc * n * x2(t)) / (kc * ka + kc * x2(t) + ka * x0(t)),
    x3'(t) = b2 * (x2(t) - x3(t)),
    y1(t) = x0(t)
)

St = StructuralIdentifiability.@ODEmodel(
    S'(t) = r * S(t) - (e + a * W(t)) * S(t) - d * W(t) * S(t) + g * R(t),
    R'(t) = rR * R(t) + (e + a * W(t)) * S(t) - dr * W(t) * R(t) - g * R(t),
    W'(t) = Dd * (T - W(t)),
    y1(t) = S(t) + R(t),
    y2(t) = T
)

cd8 = StructuralIdentifiability.@ODEmodel(
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

qwwc = StructuralIdentifiability.@ODEmodel(
    x'(t) = a * (y(t) - x(t)) + y(t) * z(t),
    y'(t) = b * (x(t) + y(t)) - x(t) * z(t),
    z'(t) = -c * z(t) - d * w(t) + x(t) * y(t),
    w'(t) = e * z(t) - f * w(t) + x(t) * y(t),
    g(t) = x(t)
)

fujita = StructuralIdentifiability.@ODEmodel(
    EGFR'(t) =
        -reaction_1_k1 * EGF_EGFR(t) + reaction_1_k2 * EGF_EGFR(t) -
        EGFR(t) * EGFR_turnover + EGFR_turnover * pro_EGFR(t),
    pAkt'(t) =
        -pAkt(t) * reaction_7_k1 - pAkt(t) * reaction_5_k1 * S6(t) +
        reaction_6_k1 * pAkt_S6(t) +
        reaction_3_k1 * pEGFR_Akt(t) +
        pAkt_S6(t) * reaction_5_k2,
    pEGFR_Akt'(t) =
        pEGFR(t) * Akt(t) * reaction_2_k1 - reaction_3_k1 * pEGFR_Akt(t) -
        pEGFR_Akt(t) * reaction_2_k2,
    S6'(t) =
        pS6(t) * reaction_8_k1 - pAkt(t) * reaction_5_k1 * S6(t) +
        pAkt_S6(t) * reaction_5_k2,
    pEGFR'(t) =
        -reaction_4_k1 * pEGFR(t) + reaction_9_k1 * EGF_EGFR(t) -
        pEGFR(t) * Akt(t) * reaction_2_k1 +
        reaction_3_k1 * pEGFR_Akt(t) +
        pEGFR_Akt(t) * reaction_2_k2,
    EGF_EGFR'(t) =
        reaction_1_k1 * EGF_EGFR(t) - reaction_9_k1 * EGF_EGFR(t) -
        reaction_1_k2 * EGF_EGFR(t),
    Akt'(t) =
        pAkt(t) * reaction_7_k1 - pEGFR(t) * Akt(t) * reaction_2_k1 +
        pEGFR_Akt(t) * reaction_2_k2,
    pAkt_S6'(t) =
        pAkt(t) * reaction_5_k1 * S6(t) - reaction_6_k1 * pAkt_S6(t) -
        pAkt_S6(t) * reaction_5_k2,
    pS6'(t) = -pS6(t) * reaction_8_k1 + reaction_6_k1 * pAkt_S6(t),
    y1(t) = pEGFR(t) * a1 + a1 * pEGFR_Akt(t),
    y2(t) = a2 * pAkt(t) + a2 * pAkt_S6(t),
    y3(t) = pS6(t) * a3
)

Bilirubin2_io = StructuralIdentifiability.@ODEmodel(
    x1'(t) =
        -(k21 + k31 + k41 + k01) * x1(t) + k12 * x2(t) + k13 * x3(t) + k14 * x4(t) + u(t),
    x2'(t) = k21 * x1(t) - k12 * x2(t),
    x3'(t) = k31 * x1(t) - k13 * x3(t),
    x4'(t) = k41 * x1(t) - k14 * x4(t),
    y1(t) = x1(t)
)

# PK1 
pk1 = StructuralIdentifiability.@ODEmodel(
    x1'(t) = u1(t)-(k1+k2)*x1(t),
    x2'(t) = k1*x1(t)-(k3+k6+k7)*x2(t)+k5*x4(t),
    x3'(t) = k2*x1(t)+k3*x2(t)-k4*x3(t),
    x4'(t) = k6*x2(t)-k5*x4(t),
    y1(t) = s2*x2(t),
    y2(t) = s3*x3(t)
)

# MAPK 5 outputs bis
MAPK_5_outputs_bis = StructuralIdentifiability.@ODEmodel(
            KS00'(t) =
                -a00 * K(t) * S00(t) +
                b00 * KS00(t) +
                gamma0100 * FS01(t) +
                gamma1000 * FS10(t) +
                gamma1100 * FS11(t),
            KS01'(t) =
                -a01 * K(t) * S01(t) + b01 * KS01(t) + c0001 * KS00(t) -
                alpha01 * F(t) * S01(t) +
                beta01 * FS01(t) +
                gamma1101 * FS11(t),
            KS10'(t) =
                -a10 * K(t) * S10(t) + b10 * KS10(t) + c0010 * KS00(t) -
                alpha10 * F(t) * S10(t) +
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
                alpha11 * F(t) * S11(t) -
                (beta11 + gamma1101 + gamma1110 + gamma1100) * FS11(t),
            S10'(t) =
                -a00 * K(t) * S00(t) + (b00 + c0001 + c0010 + c0011) * KS00(t) -
                a01 * K(t) * S01(t) + (b01 + c0111) * KS01(t) -
                a10 * K(t) * S10(t) + (b10 + c1011) * KS10(t),
            S11'(t) =
                -alpha01 * F(t) * S01(t) + (beta01 + gamma0100) * FS01(t) -
                alpha10 * F(t) * S10(t) + (beta10 + gamma1000) * FS10(t) -
                alpha11 * F(t) * S11(t) +
                (beta11 + gamma1101 + gamma1110 + gamma1100) * FS11(t),
            y0(t) = K(t),
            y1(t) = F(t),
            y2(t) = S00(t),
            y3(t) = S01(t) + S10(t),
            y4(t) = S11(t)
        )

using Nemo, Logging
Groebner = StructuralIdentifiability.Groebner
ParamPunPam = StructuralIdentifiability.ParamPunPam
Base.global_logger(ConsoleLogger(Logging.Info))

@myprof io = StructuralIdentifiability.find_ioequations(MAPK_5_outputs_bis);

@time StructuralIdentifiability.find_identifiable_functions(siwr)

@time StructuralIdentifiability.find_identifiable_functions(covid)

@myprof StructuralIdentifiability.find_identifiable_functions(covid)

fg = StructuralIdentifiability.field_generators(siwr);
id = StructuralIdentifiability.field_to_ideal(fg);

mqs = StructuralIdentifiability.IdealMQS(fg);
p = Nemo.GF(2^31 - 1)
StructuralIdentifiability.ParamPunPam.reduce_mod_p!(mqs, p)
point = rand(p, length(Nemo.gens(StructuralIdentifiability.ParamPunPam.parent_params(mqs))))
mqs_spec = StructuralIdentifiability.ParamPunPam.specialize_mod_p(mqs, point);

@time graph1, gb = Groebner.groebner_learn(mqs_spec, loglevel = -3);
@time graph2, gb_2 = Groebner.groebner_learn(mqs_spec, loglevel = 0, sweep = true);

flag1, gb1 = Groebner.groebner_apply!(graph1, mqs_spec);

@myprof Groebner.groebner_apply!(graph2, mqs_spec);

@assert flag1 && flag2 && gb1 == gb2
