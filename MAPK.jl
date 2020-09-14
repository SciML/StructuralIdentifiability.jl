include("io_equation.jl")

logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

varnames = [
    "KS00", "KS01", "KS10", "FS01", "FS10", "FS11",
    "K", "F", "S00", "S01", "S10", "S11",
    "a00", "a01", "a10", "b00", "b01", "b10", "c0001", "c0010",
    "c0011", "c0111", "c1011", "alpha01", "alpha10", "alpha11",
    "beta01", "beta10", "beta11", "gamma0100", "gamma1000", "gamma1100",
    "gamma1101", "gamma1110"
]

R, (KS00, KS01, KS10, FS01, FS10, FS11,
    K, F, S00, S01, S10, S11,
    a00, a01, a10, b00, b01, b10, c0001, c0010,
    c0011, c0111, c1011, alpha01, alpha10, alpha11,
    beta01, beta10, beta11, gamma0100, gamma1000, gamma1100,
    gamma1101, gamma1110) = PolynomialRing(QQ, varnames)

ode = ODE(Dict(
    KS00 => -a00*K*S00+b00*KS00+gamma0100*FS01+gamma1000*FS10+gamma1100*FS11,
    KS01 => -a01*K*S01+b01*KS01+c0001*KS00-alpha01*F*S01+beta01*FS01+gamma1101*FS11,
    KS10 => -a10*K*S10+b10*KS10+c0010*KS00-alpha10*F*S10+beta10*FS10+gamma1110*FS11,
    FS01 => -alpha11*F*S11+beta11*FS11+c0111*KS01+c1011*KS10+c0011*KS00,
    FS10 => a00*K*S00-(b00+c0001+c0010+c0011)*KS00,
    FS11 => a01*K*S01-(b01+c0111)*KS01,
    K => a10*K*S10-(b10+c1011)*KS10,
    F => alpha01*F*S01-(beta01+gamma0100)*FS01,
    S00 => alpha10*F*S10-(beta10+gamma1000)*FS10,
    S01 => alpha11*F*S11-(beta11+gamma1101+gamma1110+gamma1100)*FS11,
    S10 => -a00*K*S00+(b00+c0001+c0010+c0011)*KS00-a01*K*S01+(b01+c0111)*KS01-a10*K*S10+(b10+c1011)*KS10,
    S11 => -alpha01*F*S01+(beta01+gamma0100)*FS01-alpha10*F*S10+(beta10+gamma1000)*FS10-alpha11*F*S11+(beta11+gamma1101+gamma1110+gamma1100)*FS11
), [])

@time io_equations = find_ioequation(ode, [F, S00, S01, S10, S11])

print(map(length, io_equations), "\n")
