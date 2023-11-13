# A. K. Manrai and J. Gunawardena. 
# The geometry of multisite phosphorylation
# https://dx.doi.org/10.1529%2Fbiophysj.108.140632
using StructuralIdentifiability

ode = @ODEmodel(
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
    y1(t) = F(t),
    y2(t) = S00(t),
    y3(t) = S01(t),
    y4(t) = S10(t),
    y5(t) = S11(t)
)

println(assess_identifiability(ode))
