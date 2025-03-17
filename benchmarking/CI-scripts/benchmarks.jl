benchmarks = Dict(
    :LV_simple => Dict(
        :name => "Modified LV for testing",
        :ode => @ODEmodel(
            x1'(t) = (a + b) * x1(t) - c * x1(t) * x2(t),
            x2'(t) = -a * b * x2(t) + d * x1(t) * x2(t),
            y1(t) = x1(t)
        ),
    ),
    :SIWR => Dict(
        :name => "SIWR original",
        :ode => @ODEmodel(
            S'(t) = mu - bi * S(t) * I(t) - bw * S(t) * W(t) - mu * S(t) + a * R(t),
            I'(t) = bw * S(t) * W(t) + bi * S(t) * I(t) - (gam + mu) * I(t),
            W'(t) = xi * (I(t) - W(t)),
            R'(t) = gam * I(t) - (mu + a) * R(t),
            y(t) = k * I(t)
        ),
        :skip => false,
    ),
    :SIWR_simplified => Dict(
        :name => "SIWR with extra output",
        :ode => @ODEmodel(
            S'(t) = mu - bi * S(t) * I(t) - bw * S(t) * W(t) - mu * S(t) + a * R(t),
            I'(t) = bw * S(t) * W(t) + bi * S(t) * I(t) - (gam + mu) * I(t),
            W'(t) = xi * (I(t) - W(t)),
            R'(t) = gam * I(t) - (mu + a) * R(t),
            y(t) = k * I(t),
            y2(t) = S(t) + I(t) + R(t)
        ),
        :skip => false,
    ),
    :Pharm => Dict(
        :name => "Pharm",
        :ode => @ODEmodel(
            x0'(t) =
                a1 * (x1(t) - x0(t)) -
                (ka * n * x0(t)) / (kc * ka + kc * x2(t) + ka * x0(t)),
            x1'(t) = a2 * (x0(t) - x1(t)),
            x2'(t) =
                b1 * (x3(t) - x2(t)) -
                (kc * n * x2(t)) / (kc * ka + kc * x2(t) + ka * x0(t)),
            x3'(t) = b2 * (x2(t) - x3(t)),
            y1(t) = x0(t)
        ),
        :skip => false,
    ),
    :SEAIJRC => Dict(
        :name => "SEAIJRC Covid model",
        :ode => @ODEmodel(
            S'(t) = -b * S(t) * (I(t) + J(t) + q * A(t)) * Ninv,
            E'(t) = b * S(t) * (I(t) + J(t) + q * A(t)) * Ninv - k * E(t),
            A'(t) = k * (1 - r) * E(t) - g1 * A(t),
            I'(t) = k * r * E(t) - (alpha + g1) * I(t),
            J'(t) = alpha * I(t) - g2 * J(t),
            C'(t) = alpha * I(t),
            y(t) = C(t),
            y2(t) = Ninv
        ),
        :skip => false,
    ),
    :MAPK5 => Dict(
        :name => "MAPK model (5 outputs)",
        :ode => @ODEmodel(
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
            y1(t) = F(t),
            y2(t) = S00(t),
            y3(t) = S01(t),
            y4(t) = S10(t),
            y5(t) = S11(t)
        ),
    ),
    :MAPK5bis => Dict(
        :name => "MAPK model (5 outputs bis)",
        :ode => @ODEmodel(
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
        ),
        :skip => false,
    ),
    :MAPK6 => Dict(
        :name => "MAPK model (6 outputs)",
        :ode => @ODEmodel(
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
            y3(t) = S01(t),
            y4(t) = S10(t),
            y5(t) = S11(t)
        ),
    ),
    :Goodwin => Dict(
        :name => "Goodwin oscillator",
        :ode => @ODEmodel(
            x1'(t) = -b * x1(t) + 1 / (c + x4(t)),
            x2'(t) = alpha * x1(t) - beta * x2(t),
            x3'(t) = gama * x2(t) - delta * x3(t),
            x4'(t) = sigma * x4(t) * (gama * x2(t) - delta * x3(t)) / x3(t),
            y(t) = x1(t)
        ),
    ),
    :HIV => Dict(
        :name => "HIV",
        :ode => @ODEmodel(
            x'(t) = lm - d * x(t) - beta * x(t) * v(t),
            y'(t) = beta * x(t) * v(t) - a * y(t),
            v'(t) = k * y(t) - u * v(t),
            w'(t) = c * x(t) * y(t) * w(t) - c * q * y(t) * w(t) - b * w(t),
            z'(t) = c * q * y(t) * w(t) - h * z(t),
            y1(t) = w(t),
            y2(t) = z(t)
        ),
    ),
    :SIRSforced => Dict(
        :name => "SIRS forced",
        :ode => @ODEmodel(
            s'(t) = mu - mu * s(t) - b0 * (1 + b1 * x1(t)) * i(t) * s(t) + g * r(t),
            i'(t) = b0 * (1 + b1 * x1(t)) * i(t) * s(t) - (nu + mu) * i(t),
            r'(t) = nu * i(t) - (mu + g) * r(t),
            x1'(t) = -M * x2(t),
            x2'(t) = M * x1(t),
            y1(t) = i(t),
            y2(t) = r(t)
        ),
    ),
    :Akt => Dict(
        :name => "Akt pathway",
        :ode => @ODEmodel(
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
        ),
    ),
    :CD8_Tcell => Dict(
        :name => "CD8 T cell differentiation",
        :ode => @ODEmodel(
            N'(t) = -N(t) * mu_N - N(t) * P(t) * delta_NE,
            E'(t) =
                N(t) * P(t) * delta_NE - E(t)^2 * mu_EE - E(t) * delta_EL +
                E(t) * P(t) * rho_E,
            S'(t) =
                S(t) * delta_EL - S(t) * delta_LM - S(t)^2 * mu_LL - E(t) * S(t) * mu_LE,
            M'(t) = S(t) * delta_LM - mu_M * M(t),
            P'(t) =
                P(t)^2 * rho_P - P(t) * mu_P - E(t) * P(t) * mu_PE - S(t) * P(t) * mu_PL,
            y1(t) = N(t),
            y2(t) = E(t) + S(t),
            y3(t) = M(t)
        ),
    ),
    :CRN => Dict(
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
        ),
    ),
    :SLIQR => Dict(
        :name => "SLIQR",
        :ode => @ODEmodel(
            S'(t) = -b * In(t) * S(t) * Ninv - u(t) * S(t) * Ninv,
            L'(t) = b * In(t) * S(t) * Ninv - a * L(t),
            In'(t) = a * L(t) - g * In(t) + s * Q(t),
            Q'(t) = (1 - e) * g * In(t) - s * Q(t),
            y(t) = In(t) * Ninv
        ),
    ),
    :Sntg => Dict(
        :name => "Sntg",
        :ode => @ODEmodel(
            S'(t) = r * S(t) - (e + a * W(t)) * S(t) - d * W(t) * S(t) + g * R(t),
            R'(t) = rR * R(t) + (e + a * W(t)) * S(t) - dr * W(t) * R(t) - g * R(t),
            W'(t) = Dd * (T - W(t)),
            y1(t) = S(t) + R(t),
            y2(t) = T
        ),
    ),
    :QY => Dict(
        :name => "QY",
        :ode => @ODEmodel(
            P0'(t) = P1(t),
            P1'(t) = P2(t),
            P2'(t) = P3(t),
            P3'(t) = P4(t),
            P4'(t) =
                -(
                    Ks * M * siga1 * siga2 * P1(t) +
                    (
                        Ks * M * siga1 +
                        Ks * M * siga2 +
                        Ks * siga1 * siga2 +
                        siga1 * siga2 * M
                    ) * P2(t) +
                    (
                        Ks * M +
                        Ks * siga1 +
                        Ks * siga2 +
                        M * siga1 +
                        M * siga2 +
                        siga1 * siga2
                    ) * P3(t) +
                    (Ks + M + siga1 + siga2) * P4(t)
                ) -
                (
                    Mar * P5(t) +
                    beta +
                    beta_SA / (siga2 * M) * (
                        P3(t) +
                        P2(t) * (Ks + M + Mar) +
                        P1(t) * (Ks * M + Ks * Mar + M * Mar) +
                        P0(t) * Ks * M * Mar
                    ) +
                    beta_SI / M * (P2(t) + P1(t) * (Ks + Mar) + P0(t) * Ks * Mar) +
                    beta_SA * phi / ((1 - phi) * siga2 * M) * (
                        P3(t) +
                        P2(t) * (Ks + M + siga2) +
                        P1(t) * (Ks * M + Ks * siga2 + M * siga2) +
                        P0(t) * Ks * M * siga2
                    )
                ) * (
                    alpa +
                    Ks * M * siga1 * siga2 * P0(t) +
                    (
                        Ks * M * siga1 +
                        Ks * M * siga2 +
                        Ks * siga1 * siga2 +
                        siga1 * siga2 * M
                    ) * P1(t) +
                    (
                        Ks * M +
                        Ks * siga1 +
                        Ks * siga2 +
                        M * siga1 +
                        M * siga2 +
                        siga1 * siga2
                    ) * P2(t) +
                    (Ks + M + siga1 + siga2) * P3(t) +
                    P4(t)
                ),
            P5'(t) =
                -Mar * P5(t) - (
                    beta +
                    beta_SA / (siga2 * M) * (
                        P3(t) +
                        P2(t) * (Ks + M + Mar) +
                        P1(t) * (Ks * M + Ks * Mar + M * Mar) +
                        P0(t) * Ks * M * Mar
                    ) +
                    beta_SI / M * (P2(t) + P1(t) * (Ks + Mar) + P0(t) * Ks * Mar) +
                    beta_SA * phi / ((1 - phi) * siga2 * M) * (
                        P3(t) +
                        P2(t) * (Ks + M + siga2) +
                        P1(t) * (Ks * M + Ks * siga2 + M * siga2) +
                        P0(t) * Ks * M * siga2
                    )
                ),
            y(t) = P0(t)
        ),
    ),
    :LLW => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/General/LLW1987_io.jl
        :name => "LLW1987_io",
        :ode => @ODEmodel(
            x1'(t) = -p1 * x1(t) + p2 * u(t),
            x2'(t) = -p3 * x2(t) + p4 * u(t),
            x3'(t) = -(p1 + p3) * x3(t) + (p4 * x1(t) + p2 * x2(t)) * u(t),
            y1(t) = x3(t)
        ),
    ),
    :Bilirubin => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Physiology/Bilirubin2_io.jl
        :name => "Bilirubin2_io",
        :ode => @ODEmodel(
            x1'(t) =
                -(k21 + k31 + k41 + k01) * x1(t) +
                k12 * x2(t) +
                k13 * x3(t) +
                k14 * x4(t) +
                u(t),
            x2'(t) = k21 * x1(t) - k12 * x2(t),
            x3'(t) = k31 * x1(t) - k13 * x3(t),
            x4'(t) = k41 * x1(t) - k14 * x4(t),
            y1(t) = x1(t)
        ),
    ),
    # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Epidemiology/SEUIR.jl
    :SEUIR => Dict(
        :name => "SEUIR",
        :ode => @ODEmodel(
            S'(t) = -beta * (U(t) + I(t)) * (S(t) / N),
            E'(t) = beta * (U(t) + I(t)) * (S(t) / N) - E(t) * z,
            U'(t) = (z - w) * E(t) - U(t) * d,
            I'(t) = w * E(t) - I(t) * d,
            R'(t) = (U(t) + I(t)) * d,
            y1(t) = I(t)
        ),
    ),

)
