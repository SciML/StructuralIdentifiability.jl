benchmarks = Dict(
    :LV => Dict(
        :name => "Modified LV for testing",
        :ode => @ODEmodel(
            x1'(t) = (a + b) * x1(t) - c * x1(t) * x2(t),
            x2'(t) = -a * b * x2(t) + d * x1(t) * x2(t),
            y1(t) = x1(t)
        ),
    ),
    :SIWR_orig => Dict(
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
    :SIWR_multiout => Dict(
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
        :skip => false,
    ),
    :MAPK_5out => Dict(
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
    :MAPK_5out_bis => Dict(
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
    :MAPK_6out => Dict(
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
    :SIRC_forced => Dict(
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
    :CD8 => Dict(
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
    :QWWC => Dict(
        :name => "QWWC",
        :ode => @ODEmodel(
            x'(t) = a * (y(t) - x(t)) + y(t) * z(t),
            y'(t) = b * (x(t) + y(t)) - x(t) * z(t),
            z'(t) = -c * z(t) - d * w(t) + x(t) * y(t),
            w'(t) = e * z(t) - f * w(t) + x(t) * y(t),
            g(t) = x(t)
        ),
    ),
    # Gleb: 
    # The most simplified version of the result I could get by hand is
    # Ninv, b, s, g + a,  (1 - e) * g * (a - s), a * e * g
    # but I have no idea how I did this...
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
    :St => Dict(
        :name => "St",
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
    :Fujita => Dict(
        :name => "Fujita",
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
    :HIV2 => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Virology/HIV2_io.jl
        :name => "HIV2_io",
        :ode => @ODEmodel(
            x1'(t) = -b * x1(t) * x4(t) - d * x1(t) + s,
            x2'(t) = b * q1 * x1(t) * x4(t) - k1 * x2(t) - w1 * x2(t),
            x3'(t) = b * q2 * x1(t) * x4(t) + k1 * x2(t) - w2 * x3(t),
            x4'(t) = -c * x4(t) + k2 * x3(t),
            y1(t) = x1(t),
            y2(t) = x4(t)
        ),
    ),
    :Biohydrogenation => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Metabolism/Biohydrogenation_io.jl
        :name => "Biohydrogenation_io",
        :ode => @ODEmodel(
            x4'(t) = -k5 * x4(t) // (k6 + x4(t)),
            x5'(t) = k5 * x4(t) // (k6 + x4(t)) - k7 * x5(t) / (k8 + x5(t) + x6(t)),
            x6'(t) =
                k7 * x5(t) // (k8 + x5(t) + x6(t)) - k9 * x6(t) * (k10 - x6(t)) // k10,
            x7'(t) = k9 * x6(t) * (k10 - x6(t)) // k10,
            y1(t) = x4(t),
            y2(t) = x5(t)
        ),
    ),
    :Treatment => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Epidemiology/Treatment_io.jl
        :name => "Treatment_io",
        :ode => @ODEmodel(
            S'(t) = -b * S(t) * In(t) / N(t) - d * b * S(t) * Tr(t) / N(t),
            In'(t) =
                b * S(t) * In(t) / N(t) + d * b * S(t) * Tr(t) / N(t) - (a + g) * In(t),
            Tr'(t) = g * In(t) - nu * Tr(t),
            N'(t) = 0,
            y1(t) = Tr(t),
            y2(t) = N(t)
        ),
    ),
    :SEIR1 => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Epidemiology/SEIR_1_io.jl
        :name => "SEIR_1_io",
        :ode => @ODEmodel(
            S'(t) = -beta * S(t) * I(t),
            E'(t) = beta * S(t) * I(t) - v * E(t),
            I'(t) = v * E(t) - psi * I(t) - (1 - psi) * gamma * I(t),
            R'(t) = gamma * Q(t) + (1 - psi) * gamma * I(t),
            Q'(t) = -gamma * Q(t) + psi * I(t),
            y1(t) = Q(t)
        ),
    ),
    :JAK_STAT => Dict(
        # https://arxiv.org/pdf/2207.09745.pdf
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Cellular%20signalling/JAKSTAT1.jl
        :name => "JAK-STAT 1",
        :ode => @ODEmodel(
            x1'(t) = -t1 * x1(t) * 2 * u(t) - t5 * x1(t) + t6 * x2(t),
            x2'(t) = t5 * x1(t) - t6 * x2(t),
            x3'(t) = t1 * 2 * u(t) * x1(t) - t2 * x3(t) * (-x6(t) + 3),
            x4'(t) = t2 * x3(t) * (-x6(t) + 3) - t3 * x4(t),
            x5'(t) = t3 * x4(t) - t4 * x5(t),
            x6'(t) =
                -t7 * x3(t) * x6(t) / (1 + t13 * x1(t)) -
                t7 * x4(t) * x6(t) / (1 + t13 * x10(t)) + t8 * (-x6(t) + 3) * 92,
            x7'(t) = -t9 * x7(t) * (-x6(t) + 3) + t10 * (-x7(t) + 165) * 92,
            x8'(t) = t11 * (-x7(t) + 165),
            x9'(t) = -t12 * 2 * u(t) * x9(t),
            x10'(t) = x8(t) * t14 / (t15 + x8(t)) - t16 * x10(t),
            y1(t) = x1(t) + x3(t) + x4(t),
            y2(t) = t18 * (x3(t) + x4(t) + x5(t) + (1 / 3 - x9(t))),
            y3(t) = t19 * (x4(t) + x5(t)),
            y4(t) = t20 * (-x6(t) + 3),
            y5(t) = t21 * x8(t),
            y6(t) = t22 * x8(t) * t17 / t11,
            y7(t) = x10(t),
            y8(t) = -x7(t) + 165
        ),
    ),
    :SIR24 => Dict(
        # https://arxiv.org/pdf/2006.14295.pdf
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Epidemiology/SIR_24.jl
        :name => "SIR 24",
        :ode => @ODEmodel(
            A'(t) = 0,
            S'(t) = A(t) - mu + S(t) - c * phi * I(t) * S(t) / (I(t) + S(t)),
            I'(t) =
                -mu * I(t) + c * phi * I(t) * S(t) / (I(t) + S(t)) - gamma * I(t) -
                I(t) * u1(t) / (S(t) + I(t)),
            R'(t) = -mu * R(t) + gamma * I(t) + I(t) * u1(t) / (S(t) + I(t)),
            y1(t) = K * I(t)
        ),
    ),
    :SIR21 => Dict(
        # https://arxiv.org/pdf/2006.14295.pdf
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Epidemiology/SIR_21.jl
        :name => "SIR 21",
        :ode => @ODEmodel(
            S'(t) = -beta * S(t) * I(t) / N - pp * S(t) + q * C(t),
            I'(t) = beta * S(t) * I(t) / N - (r + mu) * I(t),
            R'(t) = r * I(t),
            C'(t) = pp * S(t) - q * C(t),
            D'(t) = mu * I(t),
            y1(t) = N,
            y2(t) = D(t),
            y3(t) = C(t)
        ),
    ),
    :SIR19 => Dict(
        # https://arxiv.org/pdf/2006.14295.pdf
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Epidemiology/SIR_19.jl
        :name => "SIR 19",
        :ode => @ODEmodel(
            S'(t) = -beta * S(t) * I(t) / N - pp * S(t) + q * C(t),
            I'(t) = beta * S(t) * I(t) / N - (r + mu) * I(t),
            C'(t) = pp * S(t) - q * C(t),
            D'(t) = mu * I(t),
            y1(t) = N,
            y2(t) = D(t),
            y3(t) = C(t)
        ),
    ),
    :SIR6 => Dict(
        # https://arxiv.org/pdf/2006.14295.pdf
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Epidemiology/SIR_6.jl
        :name => "SIR 6",
        :ode => @ODEmodel(
            S'(t) = -beta * I(t) * S(t) / N,
            I'(t) = beta * I(t) * S(t) / N - gamma * I(t),
            y1(t) = I(t) * K,
            y2(t) = N
        ),
    ),
    :SEIR34 => Dict(
        # https://arxiv.org/pdf/2006.14295.pdf
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Epidemiology/SEIR_34.jl
        :name => "SEIR 34",
        :ode => @ODEmodel(
            A'(t) = 0,
            S'(t) = A(t) - r * beta * S(t) * I(t) / N - mu * S(t),
            E'(t) = r * beta * S(t) * I(t) / N - epsilon * E(t) - mu * E(t),
            I'(t) = epsilon * E(t) - gamma * I(t) - mu * I(t),
            y1(t) = K * I(t),
            y2(t) = A(t),
            y3(t) = N
        ),
    ),
    :SEIR36_ref => Dict(
        # https://arxiv.org/pdf/2006.14295.pdf
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Epidemiology/SEIR_36_ref.jl
        :name => "SEIR 36 ref",
        :ode => @ODEmodel(
            S'(t) =
                -beta * S(t) * I(t) / N - q * beta_d * S(t) * Di(t) / N +
                nu * N - mu_0 * S(t),
            E'(t) =
                beta * S(t) * I(t) / N + q(t) * beta_d * S(t) * Di(t) / N - s * E(t) - phi_e * E(t) - mu_0 * E(t),
            I'(t) = s * E(t) - gamma * I(t) - mu_i * I(t) - phi * I(t) - mu_0 * I(t),
            De'(t) = phi_e * E(t) - s_d * De(t) - mu_0 * De(t),
            Di'(t) =
                phi * I(t) + s_d * De(t) - gamma_d * Di(t) - mu_d * Di(t) - mu_0 * Di(t),
            F'(t) = mu_i * I(t) + mu_d * Di(t),
            y1(t) = De(t),
            y2(t) = Di(t),
            y5(t) = F(t),
            y3(t) = N,
            y4(t) = nu,
            y6(t) = q
        ),
    ),
    :TumorHu => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Tumor/TumorHu2019.jl
        :name => "TumorHu2019",
        :ode => @ODEmodel(
            x'(t) =
                (r1 + b1 * y(t)) * x(t) * (1 - b1 * x(t)) - d1 * x(t) * z(t) / (m1 + w(t)), #pancreatic cancer cell population
            y'(t) =
                (r2 + b2 * w(t) / (k2 + w(t))) * y(t) * (1 - b2 * y(t)) - mu2 * y(t), #pancreatic stellate cell population
            z'(t) = b3 * z(t) * v(t) / ((k3 + v(t)) * (m3 + w(t))) - mu3 * z(t) + r3, #effector cells, including CD8+T cells and NK cells
            w'(t) =
                b4 * x(t) * z(t) / (k4 + x(t)) - mu4 * w(t) +
                r4 * x(t) * y(t) / (m4 + v(t)), #concentration of tumor promoting cytokines, including TGF-beta and IL-6
            v'(t) = b5 * x(t) * z(t) / (k5 + x(t)) - mu5 * v(t), #concentration of tumor suppressing cytokines, including INF-gamma and IL-2
            y1(t) = z(t)
        ),
    ),
    :TumorPillis => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Tumor/TumorPillis2007.jl
        :name => "TumorPillis2007",
        :ode => @ODEmodel(
            T'(t) =
                a * T(t) * (1 - b * T(t)) - c1 * N(t) * T(t) - D(t) * T(t) -
                KT * M(t) * T(t), #tumor cells
            L'(t) =
                -m * L(t) - q * L(t) * T(t) - ucte * L(t)^2 +
                r2 * C(t) * T(t) +
                pI * L(t) * I(t) / (gI + I(t)) +
                u1(t) - KL * M(t) * L(t), # tumor-specific effector cells, T-celss
            N'(t) =
                alpha1 - f * N(t) + g * T(t) * N(t) / (h + T(t)) - p * N(t) * T(t) -
                KN * M(t) * N(t), # non-specific effector cells, NK cells
            C'(t) = alpha2 - beta * C(t) - KC * M(t) * C(t), #circulating lymphocytes
            I'(t) =
                pt * T(t) * L(t) / (gt + T(t)) + w * L(t) * I(t) - muI * I(t) + u2(t), # IL-2, VI = u2 aplicación directa, terapia de IL2
            M'(t) = -gamma * M(t) + u1(t), #chemotherapy drug, terapia/aplicación de quimio, u1 = VM
            y1(t) = L(t),
            y2(t) = N(t),
            y3(t) = M(t)
        ),
    ),
    :SEIR2T => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Epidemiology/SEIR2T.jl
        :name => "SEIR2T",
        :ode => @ODEmodel(
            S'(t) = -b * S(t) * In(t) / N,
            E'(t) = b * S(t) * In(t) / N - nu * E(t),
            In'(t) = nu * E(t) - a * In(t),
            Cu'(t) = nu * E(t),
            y1(t) = Cu(t),
            y2(t) = N
        ),
    ),
    :SEIRT => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Epidemiology/SEIRT.jl
        :name => "SEIRT",
        :ode => @ODEmodel(
            S'(t) = -beta * I(t) * (S(t) / N),
            E'(t) = beta * I(t) * (S(t) / N) - alpha * E(t),
            I'(t) = alpha * E(t) - lambda * I(t),
            y1(t) = I(t),
            y2(t) = N
        ),
    ),
    :SEUIR => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Epidemiology/SEUIR.jl
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
    :Bruno2016 => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Gene%20expression/Bruno2016.jl
        :name => "Bruno2016",
        :ode => @ODEmodel(
            beta'(t) = -kbeta * beta(t),
            cry'(t) = -kcryOH * cry(t) - kcrybeta * cry(t),
            zea'(t) = -kzea * zea(t),
            beta10'(t) = kbeta * beta(t) + kcryOH * cry(t) - kbeta10 * beta10(t),
            OHbeta10'(t) = kcrybeta * cry(t) + kzea * zea(t) - kOHbeta10 * OHbeta10(t),
            betaio'(t) = kbeta * beta(t) + kcrybeta * cry(t) + kbeta10 * beta10(t),
            OHbetaio'(t) = kcryOH * cry(t) + kzea * zea(t) + kOHbeta10 * OHbeta10(t),
            y1(t) = beta(t),
            y2(t) = beta10(t)
        ),
    ),
    :Transfection => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Gene%20expression/Transfection_4State.jl
        :name => "Transfection_4State",
        :ode => @ODEmodel(
            mRNA'(t) = -d1 * mRNA(t) - d2 * mRNA(t) * enz(t),
            GFP'(t) = kTL * mRNA(t) - b * GFP(t),
            enz'(t) = d3 * mRNAenz(t) - d2 * mRNA(t) * enz(t),
            mRNAenz'(t) = -d3 * mRNAenz(t) + d2 * mRNA(t) * enz(t),
            y1(t) = GFP(t)
        ),
    ),
    :p53 => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Gene%20expression/p53.jl
        :name => "p53",
        :ode => @ODEmodel(
            x1'(t) =
                (p1 * x4(t)) - (p3 * x1(t)) -
                p4 * ((x1(t)^2 / (p5 + x1(t))) * (1 + (p6 * u1(t) / (p7 + u1(t))))),
            x2'(t) =
                p8 - (p9 * x2(t)) -
                p10 * (
                    (x1(t) * x2(t) / (p11 + x2(t))) * (1 + (p12 * u1(t) / (p13 + u1(t))))
                ),
            x3'(t) =
                p14 - (p15 * x3(t)) -
                p16 * x1(t) * x3(t) * (1 - p18 * u1(t)) / (p17 + x3(t)),
            x4'(t) =
                p20 - p21 * (1 - p24) * (1 - p25) / ((p22^4) + 1) - (p20 * x4(t)) +
                (p21 * (x3(t)^4)) +
                (1 + p23 * u1(t)) * (1 - p24 * x1(t)) * (1 - p25 * x2(t)) /
                (p22^4 + x3(t)^4),
            y1(t) = x1(t),
            y2(t) = x2(t),
            y3(t) = x3(t),
            y4(t) = x4(t)
        ),
    ),
    :Crauste => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/General/Crauste_SI.jl
        :name => "Crauste_SI",
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
    :HDNL => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/General/HighDimNonLin.jl
        :name => "HighDimNonLin",
        :ode => @ODEmodel(
            x1'(t) = -vm * x1(t) / (km + x1(t)) - p1 * x1(t) + u(t),
            x2'(t) = p1 * x1(t) - p2 * x2(t),
            x3'(t) = p2 * x2(t) - p3 * x3(t),
            x4'(t) = p3 * x3(t) - p4 * x4(t),
            x5'(t) = p4 * x4(t) - p5 * x5(t),
            x6'(t) = p5 * x5(t) - p6 * x6(t),
            x7'(t) = p6 * x6(t) - p7 * x7(t),
            x8'(t) = p7 * x7(t) - p8 * x8(t),
            x9'(t) = p8 * x8(t) - p9 * x9(t),
            x10'(t) = p9 * x9(t) - p10 * x10(t),
            x11'(t) = p10 * x10(t) - p11 * x11(t),
            x12'(t) = p11 * x11(t) - p12 * x12(t),
            x13'(t) = p12 * x12(t) - p13 * x13(t),
            x14'(t) = p13 * x13(t) - p14 * x14(t),
            x15'(t) = p14 * x14(t) - p15 * x15(t),
            x16'(t) = p15 * x15(t) - p16 * x16(t),
            x17'(t) = p16 * x16(t) - p17 * x17(t),
            x18'(t) = p17 * x17(t) - p18 * x18(t),
            x19'(t) = p18 * x18(t) - p19 * x19(t),
            x20'(t) = p19 * x19(t) - p20 * x20(t),
            y1(t) = x1(t),
            y2(t) = x2(t),
            y3(t) = x3(t),
            y4(t) = x4(t),
            y5(t) = x5(t),
            y6(t) = x6(t),
            y7(t) = x7(t),
            y8(t) = x8(t),
            y9(t) = x9(t),
            y10(t) = x10(t),
            y11(t) = x11(t),
            y12(t) = x12(t),
            y13(t) = x13(t),
            y14(t) = x14(t),
            y15(t) = x15(t),
            y16(t) = x16(t),
            y17(t) = x17(t),
            y18(t) = x18(t),
            y19(t) = x19(t),
            y20(t) = x20(t)
        ),
    ),
    :KD1999 => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/General/KD1999.jl
        :name => "KD1999",
        :ode => @ODEmodel(
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
        ),
    ),
    :CGV1990 => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Immunology/CGV1990.jl
        :name => "CGV1990",
        :ode => @ODEmodel(
            q1'(t) = k4 * q3(t) - (k3 + k7) * q1(t) + u(t),
            q3'(t) =
                k3 * q1(t) - k4 * q3(t) - k5 * q3(t) * (R * V3 - q35(t)) + k6 * q35(t) - k5 * q3(t) * (5 * V36 / V3) * (S * V36 - q36(t)) +
                k6 * q36(t),
            q35'(t) = k5 * q3(t) * (R * V3 - q35(t)) - k6 * q35(t),
            q36'(t) = k5 * q3(t) * (5 * V36 / V3) * (S * V36 - q36(t)) - k6 * q36(t),
            q7'(t) = k7 * q1(t),
            y1(t) = q7(t)
        ),
    ),
    :LeukaemiaLeon2021 => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Immunology/LeukaemiaLeon2021.jl
        :name => "LeukaemiaLeon2021",
        :ode => @ODEmodel(
            C'(t) = rhoc * (L(t) + B(t)) * C(t) + rhob * I(t) * C(t) - C(t) / taoc, #number of cells of CAR T cells
            L'(t) = rhol * L(t) - alpha * L(t) * C(t), # number of cells of leukaemic cells
            B'(t) = I(t) / taoi - alpha * B(t) * C(t) - B(t) / taob, # number of cells of mature healthy B cells
            P'(t) =
                rhop * (2 * ap * (1 / (1 + ks * (P(t) + I(t)))) - 1) * P(t) - P(t) / taop, # number of cells of CD19- haematopoietic stem cells (HSCs)
            I'(t) =
                rhoi * (2 * ai * (1 / (1 + ks * (P(t) + I(t)))) - 1) * I(t) - I(t) / taoi + P(t) / taop - alpha * beta * I(t) * C(t), # number of cells of CD19+ B cell progenitors
            y1(t) = C(t),
            y2(t) = L(t),
            y3(t) = P(t) + B(t) + I(t) #células B totales 
        ),
    ),
    :Lipolysis => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Metabolism/Ruminal%20lipolysis.jl
        :name => "Ruminal lipolysis",
        :ode => @ODEmodel(
            x1'(t) = -x1(t) * x5(t) / (k2 + x1(t)),
            x2'(t) = 2 * x1(t) * x5(t) / ((k2 + x1(t)) * 3) - k4 * x2(t),
            x3'(t) = k4 * (x2(t)) / 2 - k4 * x3(t),
            x4'(t) = x1(t) * x5(t) / (3 * (k2 + x1(t))) + k4 * (x2(t)) / 2 + k4 * x3(t),
            x5'(t) = -k3 * x5(t),
            y1(t) = x1(t),
            y2(t) = x2(t) + x3(t),
            y3(t) = x4(t)
        ),
    ),
    :cLV1 => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Microbial/cLV1.jl
        :name => "cLV1 (2o)",
        :ode => @ODEmodel(
            pi1'(t) =
                pi1(t) * (
                    (g1 + A11 * pi1(t) + A12 * pi2(t) + A13 * pi3(t) + B11 * u1(t)) -
                    pi1(t) * (
                        g1 + A11 * pi1(t) + A12 * pi2(t) + A13 * pi3(t) + B11 * u1(t)
                    ) +
                    pi2(t) *
                    (g2 + A21 * pi1(t) + A22 * pi2(t) + A23 * pi3(t) + B21 * u1(t))
                ),
            pi2'(t) =
                pi2(t) * (
                    (g2 + A21 * pi1(t) + A22 * pi2(t) + A23 * pi3(t) + B21 * u1(t)) -
                    pi1(t) * (
                        g1 + A11 * pi1(t) + A12 * pi2(t) + A13 * pi3(t) + B11 * u1(t)
                    ) +
                    pi2(t) *
                    (g2 + A21 * pi1(t) + A22 * pi2(t) + A23 * pi3(t) + B21 * u1(t))
                ),
            pi3'(t) =
                pi3(t) * (
                    (g3 + A31 * pi1(t) + A32 * pi2(t) + A33 * pi3(t) + B31 * u1(t)) -
                    pi1(t) * (
                        g1 + A11 * pi1(t) + A12 * pi2(t) + A13 * pi3(t) + B11 * u1(t)
                    ) +
                    pi2(t) *
                    (g2 + A21 * pi1(t) + A22 * pi2(t) + A23 * pi3(t) + B21 * u1(t))
                ),
            y1(t) = pi1(t),
            y2(t) = pi2(t)
        ),
    ),
    :cLV2 => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Microbial/cLV1.jl
        :name => "cLV1 (1o)",
        :ode => @ODEmodel(
            pi1'(t) =
                pi1(t) * (
                    (g1 + A11 * pi1(t) + A12 * pi2(t) + A13 * pi3(t) + B11 * u1(t)) -
                    pi1(t) * (
                        g1 + A11 * pi1(t) + A12 * pi2(t) + A13 * pi3(t) + B11 * u1(t)
                    ) +
                    pi2(t) *
                    (g2 + A21 * pi1(t) + A22 * pi2(t) + A23 * pi3(t) + B21 * u1(t))
                ),
            pi2'(t) =
                pi2(t) * (
                    (g2 + A21 * pi1(t) + A22 * pi2(t) + A23 * pi3(t) + B21 * u1(t)) -
                    pi1(t) * (
                        g1 + A11 * pi1(t) + A12 * pi2(t) + A13 * pi3(t) + B11 * u1(t)
                    ) +
                    pi2(t) *
                    (g2 + A21 * pi1(t) + A22 * pi2(t) + A23 * pi3(t) + B21 * u1(t))
                ),
            pi3'(t) =
                pi3(t) * (
                    (g3 + A31 * pi1(t) + A32 * pi2(t) + A33 * pi3(t) + B31 * u1(t)) -
                    pi1(t) * (
                        g1 + A11 * pi1(t) + A12 * pi2(t) + A13 * pi3(t) + B11 * u1(t)
                    ) +
                    pi2(t) *
                    (g2 + A21 * pi1(t) + A22 * pi2(t) + A23 * pi3(t) + B21 * u1(t))
                ),
            y1(t) = pi1(t)
        ),
    ),
    :genLV => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Microbial/generalizedLoktaVolterra.jl
        :name => "generalizedLoktaVolterra (1o)",
        :ode => @ODEmodel(
            x1'(t) = r1 * x1(t) + beta11 * x1(t)^2 + beta12 * x1(t) * x2(t),
            x2'(t) = r2 * x2(t) + beta21 * x1(t) * x2(t) + beta22 * x2(t)^2,
            y1(t) = x1(t)
        ),
    ),
    :Pivastatin => Dict(
        # https://github.com/Xabo-RB/Local-Global-Models/blob/main/Models/Pharmacokinetics/Pivastatin.jl
        :name => "Pivastatin",
        :ode => @ODEmodel(
            x1'(t) = k3 * x3(t) - r3 * x1(t) - k1 * x1(t) * (T0 - x2(t)) + r1 * x2(t),
            x2'(t) = k1 * x1(t) * (T0 - x2(t)) - (r1 + k2) * x2(t),
            x3'(t) = r3 * x1(t) - (k3 + k4) * x3(t) + k2 * x2(t),
            y1(t) = k * (x2(t) + x3(t))
        ),
    ),
    # Next two models are from Appendix A from https://arxiv.org/abs/2412.05283
    :Lincomp1 => Dict(
        :name => "Linear_compartment_hard_1",
        :ode => linear_compartment_model([[2], [3], [4], [5], [1]], inputs = [1], outputs = [5], leaks = [2, 3]),
    ),
    :Lincomp2 => Dict(
        :name => "Linear_compartment_hard_2",
        :ode => linear_compartment_model([[2], [3], [4], [5], [1]], inputs = [1], outputs = [5], leaks = [2, 4]),
    ),
    # Equations (8)-(13) from https://www.sciencedirect.com/science/article/pii/S0022519320303945?via%3Dihub
    :Covid1 => Dict(
        :name => "Covid model (Gevertz et al)",
        :ode => @ODEmodel(
            Sd'(t) = -es * ba * (An(t) + ea * Ad(t)) * Sd(t) - h1 * Sd(t) + h2 * Sn(t) - es * bi * Sd(t) * I(t),
            Sn'(t) = -bi * Sn(t) * I(t) - ba * (An(t) + ea * Ad(t)) * Sn(t) + h1 * Sd(t) - h2 * Sn(t),
            Ad'(t) = es * bi * Sd(t) * I(t) + es * ba * (An(t) + ea * Ad(t)) * Sd(t) + h2 * An(t) - gai * Ad(t) - h1 * Ad(t),
            An'(t) = bi * Sn(t) * I(t) + ba * (An(t) + ea * Ad(t)) * Sn(t) + h1 * Ad(t) - gai * An(t) - h2 * An(t),
            I'(t) = f * gai * (Ad(t) + An(t)) - delta * I(t) - gir * I(t),
            y1(t) = Sd(t),
            y2(t) = I(t),
        ),
    ),
   # Equations (2.17)-(2.22) from https://thesis.unipd.it/bitstream/20.500.12608/15925/1/tesi.pdf
    :Covid2 => Dict(
        :name => "Covid model (Gallina)",
        :ode => @ODEmodel(
            S'(t) = b * N - S(t) * (l * I(t) + l * ea * eq * Q(t) + l * ea * A(t) + l * ej * J(t) + d1),
            I'(t) = k1 * A(t) - (g1 + m2 + d2) * I(t),
            A'(t) = S(t) * (l * I(t) + l * ea * eq * Q(t) + l * ea * A(t) + l * ej * J(t)) - (k1 + m1 + d4) * A(t),
            Q'(t) = m1 * A(t) - (k2 + d5) * Q(t),
            J'(t) = k2 * Q(t) + m2 * I(t) - (g2 + d6) * J(t),
            y1(t) = Q(t),
            y2(t) = J(t),
        ),
    ),
    # Slide 35 from https://indico.ictp.it/event/7960/session/3/contribution/19/material/slides/0.pdf
    :Covid3 => Dict(
        :name => "Covid model (Chitnis)",
        :ode => @ODEmodel(
            S'(t) = L - r * b * S(t) * I(t) / N - m * S(t),
            E'(t) = b * S(t) * I(t) / N - e * E(t) - m * E(t),
            I'(t) = e * E(t) - g * I(t) - m * I(t),
            y(t) = K * I(t),
        )
    ),
)

# the NFkB example

ode = @ODEmodel(
    x1'(t) = k_prod - k_deg * x1(t) - k1 * x1(t) * u(t),
    x2'(t) =
        -k3 * x2(t) - k_deg * x2(t) - a2 * x2(t) * x10(t) + t1 * x4(t) -
        a3 * x2(t) * x13(t) +
        t2 * x5(t) +
        (k1 * x1(t) - k2 * x2(t) * x8(t)) * u(t),
    x3'(t) = k3 * x2(t) - k_deg * x3(t) + k2 * x2(t) * x8(t) * u(t),
    x4'(t) = a2 * x2(t) * x10(t) - t1 * x4(t),
    x5'(t) = a3 * x2(t) * x13(t) - t2 * x5(t),
    x6'(t) = c6a * x13(t) - a1 * x6(t) * x10(t) + t2 * x5(t) - i1 * x6(t),
    x7'(t) = i1 * kv * x6(t) - a1 * x11(t) * x7(t),
    x8'(t) = c4 * x9(t) - c5 * x8(t),
    x9'(t) = c2 - c1 * x7(t) - c3 * x9(t),
    x10'(t) =
        -a2 * x2(t) * x10(t) - a1 * x10(t) * x6(t) + c4a * x12(t) - c5a * x10(t) -
        i1a * x10(t) + e1a * x11(t),
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

ode = set_parameter_values(
    ode,
    Dict(
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
        c3c => QQ(4, 10^(4)),
    ),
)

benchmarks[:NFkB] = Dict(:name => "NFkB", :ode => ode, :skip => false)
