benchmarks = [
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
            y(t) = C(t),
            y2(t) = Ninv(t)
        ),
        :skip => true
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
        :name => "Goodwin oscillator",
        :ode => @ODEmodel(
            x1'(t) = -b * x1(t) + 1 / (c + x4(t)),
            x2'(t) = alpha * x1(t) - beta * x2(t),
            x3'(t) = gama * x2(t) - delta * x3(t),
            x4'(t) = sigma * x4(t) * (gama * x2(t) - delta * x3(t)) / x3(t),
            y(t) = x1(t)
        )
    )
]
