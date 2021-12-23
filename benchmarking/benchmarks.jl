benchmarks = [
    Dict(
        :name => "Modified LV for testing",
        :ode => @ODEmodel(
            x1'(t) = (a + b) * x1(t) - c * x1(t) * x2(t),
            x2'(t) = -a * b * x2(t) + d * x1(t) * x2(t),
            y1(t) = x1(t)
        )
    ),
    Dict(
        :name => "SIWR original",
        :ode => @ODEmodel(
            S'(t) = mu - bi * S(t) * I(t) - bw * S(t) * W(t) - mu * S(t) + a * R(t),
            I'(t) = bw * S(t) * W(t) + bi * S(t) * I(t) - (gam + mu) * I(t),
            W'(t) = xi * (I(t) - W(t)),
            R'(t) = gam * I(t) - (mu + a) * R(t),
            y(t) = k * I(t)
        ),
        :skip => false
    ),
    Dict(
        :name => "SIWR with extra output",
        :ode => @ODEmodel(
            S'(t) = mu - bi * S(t) * I(t) - bw * S(t) * W(t) - mu * S(t) + a * R(t),
            I'(t) = bw * S(t) * W(t) + bi * S(t) * I(t) - (gam + mu) * I(t),
            W'(t) = xi * (I(t) - W(t)),
            R'(t) = gam * I(t) - (mu + a) * R(t),
            y(t) = k * I(t),
            y2(t) = S(t) + I(t) + R(t)
        ),
        :skip => false
    ),
    Dict(
        :name => "Pharm",
        :ode => @ODEmodel(
            x0'(t) = a1 * (x1(t) - x0(t)) - (ka * n * x0(t)) / (kc * ka + kc * x2(t) + ka * x0(t)),
            x1'(t) = a2 * (x0(t) - x1(t)),
            x2'(t) = b1 * (x3(t) - x2(t)) - (kc * n * x2(t)) / (kc * ka + kc * x2(t) + ka * x0(t)),
            x3'(t) = b2 * (x2(t) - x3(t)),
            y1(t) = x0(t)
        ),
        :skip => false
    ),
    Dict(
        :name => "SEAIJRC Covid model",
        :ode => @ODEmodel(
            S'(t) = -b * S(t) * (I(t) + J(t) + q * A(t)) * Ninv(t),
            E'(t) = b * S(t) * (I(t) + J(t) + q * A(t)) * Ninv(t) - k * E(t),
            A'(t) = k * (1 - r) * E(t) - g1 * A(t),
            I'(t) = k * r * E(t) - (alpha + g1) * I(t),
            J'(t) = alpha * I(t) - g2 * J(t),
            C'(t) = alpha * I(t),
            Ninv'(t) = 0,
            y(t) = C(t),
            y2(t) = Ninv(t)
        ),
        :skip => false
    ),
    Dict(
        :name => "MAPK model (5 outputs)",
        :ode => @ODEmodel(
            KS00'(t) = -a00 * K(t) * S00(t) + b00 * KS00(t) + gamma0100 * FS01(t) + gamma1000 * FS10(t) + gamma1100 * FS11(t),
            KS01'(t) = -a01 * K(t) * S01(t) + b01 * KS01(t) + c0001 * KS00(t) - alpha01 * F(t) * S01(t) + beta01 * FS01(t) + gamma1101 * FS11(t),
            KS10'(t) = -a10 * K(t) * S10(t) + b10 * KS10(t) + c0010 * KS00(t) - alpha10 * F(t) * S10(t) + beta10 * FS10(t) + gamma1110 * FS11(t),
            FS01'(t) = -alpha11 * F(t) * S11(t) + beta11 * FS11(t) + c0111 * KS01(t) + c1011 * KS10(t) + c0011 * KS00(t),
            FS10'(t) = a00 * K(t) * S00(t) - (b00 + c0001 + c0010 + c0011) * KS00(t),
            FS11'(t) = a01 * K(t) * S01(t) - (b01 + c0111) * KS01(t),
            K'(t) = a10 * K(t) * S10(t) - (b10 + c1011) * KS10(t),
            F'(t) = alpha01 * F(t) * S01(t) - (beta01 + gamma0100) * FS01(t),
            S00'(t) = alpha10 * F(t) * S10(t) - (beta10 + gamma1000) * FS10(t),
            S01'(t) = alpha11 * F(t) * S11(t) - (beta11 + gamma1101 + gamma1110 + gamma1100) * FS11(t),
            S10'(t) = -a00 * K(t) * S00(t) + (b00 + c0001 + c0010 + c0011) * KS00(t) - a01 * K(t) * S01(t) + (b01 + c0111) * KS01(t) - a10 * K(t) * S10(t) + (b10 + c1011) * KS10(t),
            S11'(t) = -alpha01 * F(t) * S01(t) + (beta01 + gamma0100) * FS01(t) - alpha10 * F(t) * S10(t) + (beta10 + gamma1000) * FS10(t) - alpha11 * F(t) * S11(t) + (beta11 + gamma1101 + gamma1110 + gamma1100) * FS11(t),
            y1(t) = F(t),
            y2(t) = S00(t),
            y3(t) = S01(t),
            y4(t) = S10(t),
            y5(t) = S11(t)
        ) 
    ),
    Dict(
        :name => "MAPK model (5 outputs bis)",
        :ode => @ODEmodel(
            KS00'(t) = -a00 * K(t) * S00(t) + b00 * KS00(t) + gamma0100 * FS01(t) + gamma1000 * FS10(t) + gamma1100 * FS11(t),
            KS01'(t) = -a01 * K(t) * S01(t) + b01 * KS01(t) + c0001 * KS00(t) - alpha01 * F(t) * S01(t) + beta01 * FS01(t) + gamma1101 * FS11(t),
            KS10'(t) = -a10 * K(t) * S10(t) + b10 * KS10(t) + c0010 * KS00(t) - alpha10 * F(t) * S10(t) + beta10 * FS10(t) + gamma1110 * FS11(t),
            FS01'(t) = -alpha11 * F(t) * S11(t) + beta11 * FS11(t) + c0111 * KS01(t) + c1011 * KS10(t) + c0011 * KS00(t),
            FS10'(t) = a00 * K(t) * S00(t) - (b00 + c0001 + c0010 + c0011) * KS00(t),
            FS11'(t) = a01 * K(t) * S01(t) - (b01 + c0111) * KS01(t),
            K'(t) = a10 * K(t) * S10(t) - (b10 + c1011) * KS10(t),
            F'(t) = alpha01 * F(t) * S01(t) - (beta01 + gamma0100) * FS01(t),
            S00'(t) = alpha10 * F(t) * S10(t) - (beta10 + gamma1000) * FS10(t),
            S01'(t) = alpha11 * F(t) * S11(t) - (beta11 + gamma1101 + gamma1110 + gamma1100) * FS11(t),
            S10'(t) = -a00 * K(t) * S00(t) + (b00 + c0001 + c0010 + c0011) * KS00(t) - a01 * K(t) * S01(t) + (b01 + c0111) * KS01(t) - a10 * K(t) * S10(t) + (b10 + c1011) * KS10(t),
            S11'(t) = -alpha01 * F(t) * S01(t) + (beta01 + gamma0100) * FS01(t) - alpha10 * F(t) * S10(t) + (beta10 + gamma1000) * FS10(t) - alpha11 * F(t) * S11(t) + (beta11 + gamma1101 + gamma1110 + gamma1100) * FS11(t),
            y0(t) = K(t),
            y1(t) = F(t),
            y2(t) = S00(t),
            y3(t) = S01(t) + S10(t),
            y4(t) = S11(t)
        ),
        :skip => false 
    ),
    Dict(
        :name => "MAPK model (6 outputs)",
        :ode => @ODEmodel(
            KS00'(t) = -a00 * K(t) * S00(t) + b00 * KS00(t) + gamma0100 * FS01(t) + gamma1000 * FS10(t) + gamma1100 * FS11(t),
            KS01'(t) = -a01 * K(t) * S01(t) + b01 * KS01(t) + c0001 * KS00(t) - alpha01 * F(t) * S01(t) + beta01 * FS01(t) + gamma1101 * FS11(t),
            KS10'(t) = -a10 * K(t) * S10(t) + b10 * KS10(t) + c0010 * KS00(t) - alpha10 * F(t) * S10(t) + beta10 * FS10(t) + gamma1110 * FS11(t),
            FS01'(t) = -alpha11 * F(t) * S11(t) + beta11 * FS11(t) + c0111 * KS01(t) + c1011 * KS10(t) + c0011 * KS00(t),
            FS10'(t) = a00 * K(t) * S00(t) - (b00 + c0001 + c0010 + c0011) * KS00(t),
            FS11'(t) = a01 * K(t) * S01(t) - (b01 + c0111) * KS01(t),
            K'(t) = a10 * K(t) * S10(t) - (b10 + c1011) * KS10(t),
            F'(t) = alpha01 * F(t) * S01(t) - (beta01 + gamma0100) * FS01(t),
            S00'(t) = alpha10 * F(t) * S10(t) - (beta10 + gamma1000) * FS10(t),
            S01'(t) = alpha11 * F(t) * S11(t) - (beta11 + gamma1101 + gamma1110 + gamma1100) * FS11(t),
            S10'(t) = -a00 * K(t) * S00(t) + (b00 + c0001 + c0010 + c0011) * KS00(t) - a01 * K(t) * S01(t) + (b01 + c0111) * KS01(t) - a10 * K(t) * S10(t) + (b10 + c1011) * KS10(t),
            S11'(t) = -alpha01 * F(t) * S01(t) + (beta01 + gamma0100) * FS01(t) - alpha10 * F(t) * S10(t) + (beta10 + gamma1000) * FS10(t) - alpha11 * F(t) * S11(t) + (beta11 + gamma1101 + gamma1110 + gamma1100) * FS11(t),
            y0(t) = K(t),
            y1(t) = F(t),
            y2(t) = S00(t),
            y3(t) = S01(t),
            y4(t) = S10(t),
            y5(t) = S11(t)
        ) 
    ),
    Dict(
        :name => "Goodwin oscillator",
        :ode => @ODEmodel(
            x1'(t) = -b * x1(t) + 1 / (c + x4(t)),
            x2'(t) = alpha * x1(t) - beta * x2(t),
            x3'(t) = gama * x2(t) - delta * x3(t),
            x4'(t) = sigma * x4(t) * (gama * x2(t) - delta * x3(t)) / x3(t),
            y(t) = x1(t)
        )
    ),
    Dict(
        :name => "HIV",
        :ode => @ODEmodel(
            x'(t) = lm - d * x(t) - beta * x(t) * v(t),
            y'(t) = beta * x(t) * v(t) - a * y(t),
            v'(t) = k * y(t) - u * v(t),
            w'(t) = c * x(t) * y(t) * w(t) - c * q * y(t) * w(t) - b * w(t),
            z'(t) = c * q * y(t) * w(t) - h * z(t),
            y1(t) = w(t),
            y2(t) = z(t)
        )
    ),
    Dict(
        :name => "SIRS forced",
        :ode => @ODEmodel(
            s'(t) = mu - mu * s(t) - b0 * (1 + b1 * x1(t)) * i(t) * s(t) + g * r(t),
            i'(t) = b0 * (1 + b1 * x1(t)) * i(t) * s(t) - (nu + mu) * i(t),
            r'(t) = nu * i(t) - (mu + g) * r(t),
            x1'(t) = -M * x2(t),
            x2'(t) = M * x1(t),
            y1(t) = i(t),
            y2(t) = r(t)
        )
    ),
    Dict(
        :name => "Akt pathway",
        :ode => @ODEmodel(
            EGFR'(t) = EGFR_turnover * pro_EGFR(t) + EGF_EGFR(t) * reaction_1_k2 - EGFR(t) * EGFR_turnover - EGF_EGFR(t) * reaction_1_k1,
            pEGFR'(t) = EGF_EGFR(t) * reaction_9_k1 - pEGFR(t) * reaction_4_k1 + pEGFR_Akt(t) * reaction_2_k2 + pEGFR_Akt(t) * reaction_3_k1 - Akt(t) * pEGFR(t) * reaction_2_k1,
            pEGFR_Akt'(t) = Akt(t) * pEGFR(t) * reaction_2_k1 - pEGFR_Akt(t) * reaction_3_k1 - pEGFR_Akt(t) * reaction_2_k2,
            Akt'(t) = pAkt(t) * reaction_7_k1 + pEGFR_Akt(t) * reaction_2_k2 - Akt(t) * pEGFR(t) * reaction_2_k1,
            pAkt'(t) = pAkt_S6(t) * reaction_5_k2 - pAkt(t) * reaction_7_k1 + pAkt_S6(t) * reaction_6_k1 + pEGFR_Akt(t) * reaction_3_k1 - S6(t) * pAkt(t) * reaction_5_k1,
            S6'(t) = pAkt_S6(t) * reaction_5_k2 + pS6(t) * reaction_8_k1 - S6(t) * pAkt(t) * reaction_5_k1,
            pAkt_S6'(t) = S6(t) * pAkt(t) * reaction_5_k1 - pAkt_S6(t) * reaction_6_k1 - pAkt_S6(t) * reaction_5_k2,
            pS6'(t) = pAkt_S6(t) * reaction_6_k1 - pS6(t) * reaction_8_k1,
            EGF_EGFR'(t) = EGF_EGFR(t) * reaction_1_k1 - EGF_EGFR(t) * reaction_9_k1 - EGF_EGFR(t) * reaction_1_k2,
            y1(t) = a1 * (pEGFR(t) + pEGFR_Akt(t)),
            y2(t) = a2 * (pAkt(t) + pAkt_S6(t)),
            y3(t) = a3 * pS6(t)
        )
    ),
    Dict(
        :name => "CD8 T cell differentiation",
        :ode => @ODEmodel(
            N'(t) = - N(t) * mu_N - N(t) * P(t) * delta_NE,
            E'(t) = N(t) * P(t) * delta_NE - E(t)^2 * mu_EE - E(t) * delta_EL + E(t) * P(t) * rho_E,
            S'(t) = S(t) * delta_EL - S(t) * delta_LM - S(t)^2 * mu_LL - E(t) * S(t) * mu_LE,
            M'(t) = S(t) * delta_LM - mu_M * M(t),
            P'(t) = P(t)^2 * rho_P - P(t) * mu_P - E(t) * P(t) * mu_PE - S(t) * P(t) * mu_PL,
            y1(t) = N(t),
            y2(t) = E(t) + S(t),
            y3(t) = M(t)
        )
    ),
    Dict(
        :name => "Chemical reaction network",
        :ode => @ODEmodel(    
            x1'(t) = -k1 * x1(t) * x2(t) + k2 * x4(t) + k4 * x6(t),
            x2'(t) = -k1 * x1(t) * x2(t) + k2 * x4(t) + k3 * x4(t),
            x3'(t) = k3 * x4(t) + k5 * x6(t) - k6 * x3(t) * x5(t),
            x4'(t) = k1 * x1(t) * x2(t) - k2 * x4(t) - k3 * x4(t),
            x5'(t) = k4 * x6(t) + k5 * x6(t) - k6 * x3(t) * x5(t),
            x6'(t) = -k4 * x6(t) - k5 * x6(t) + k6 * x3(t) * x5(t),
            y1(t) = x3(t),
            y2(t) = x2(t)
        )
    )
]

# the NFkB example

ode = @ODEmodel(
    x1'(t) = k_prod  - k_deg * x1(t) - k1 * x1(t) * u(t),
    x2'(t) = -k3 * x2(t) - k_deg * x2(t) - a2 * x2(t) * x10(t) + t1 * x4(t) - a3 * x2(t) * x13(t) + t2 * x5(t) + (k1 * x1(t) - k2 * x2(t) * x8(t)) * u(t),
    x3'(t) = k3 * x2(t) - k_deg * x3(t) + k2 * x2(t) * x8(t) * u(t),
    x4'(t) = a2 * x2(t) * x10(t) - t1 * x4(t),
    x5'(t) = a3 *  x2(t) * x13(t) - t2 *  x5(t),
    x6'(t) = c6a * x13(t) - a1 * x6(t) * x10(t) + t2 * x5(t) - i1 * x6(t),
    x7'(t) = i1 * kv * x6(t) - a1 * x11(t) *  x7(t),
    x8'(t) = c4 * x9(t) - c5 * x8(t),
    x9'(t) = c2 - c1 * x7(t) - c3 * x9(t),
    x10'(t) = -a2 * x2(t) * x10(t) - a1 * x10(t) * x6(t) + c4a * x12(t) - c5a * x10(t) - i1a * x10(t) + e1a * x11(t),
    x11'(t) = -a1 * x11(t) * x7(t) + i1a * kv * x10(t) - e1a * kv * x11(t),
    x12'(t) = c2a + c1a * x7(t) - c3a * x12(t),
    x13'(t) = a1 * x10(t) * x6(t) - c6a * x13(t) - a3 * x2(t) * x13(t) + e2a * x14(t),
    x14'(t) = a1 * x11(t) * x7(t) - e2a * kv * x14(t),
    x15'(t) = c2c + c1c * x7(t) - c3c * x15(t),
    y1(t) = x7(t),
    y2(t) = x10(t) + x13(t),
    y3(t) = x9(t),
    y4(t) = x1(t) + x2(t) + x3(t),
    y5(t) = x2(t),
    y6(t) = x12(t)
)

QQ = StructuralIdentifiability.Nemo.QQ

ode = set_parameter_values(ode, Dict(
    a1 => QQ(1, 2),
    a2 => QQ(1, 5),
    a3 => QQ(1),
    c1a => QQ(5, 10^(7)),
    c2a => QQ(0),
    c5a => QQ(1, 10^(4)),
    c6a => QQ(2, 10^(5)),
    c1 => QQ(5, 10^(7)),
    c2 => QQ(0),
    c3 => QQ(4, 10^(4)),
    c4 => QQ(1, 2),
    kv => QQ(5),
    e1a => QQ(5, 10^(4)),
    c1c => QQ(5, 10^(7)),
    c2c => QQ(0),
    c3c => QQ(4, 10^(4))                             
))

push!(
    benchmarks, 
    Dict(
        :name => "NFkB",
        :ode => ode,
        :skip => true
    )
)
