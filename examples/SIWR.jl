include("../io_equation.jl")

# SIWR Cholera model
logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(
    S'(t) = mu - bi * S(t) * I(t) - bw * S(t) * W(t) - mu * S(t) + a * R(t),
    I'(t) = bw * S(t) * W(t) + bi * S(t) * I(t) - (gam + mu) * I(t),
    W'(t) = xi * (I(t) - W(t)),
    R'(t) = gam * I(t) - (mu + a) * R(t),
    [k]
)

@time io_equation = collect(values(find_ioequations(ode, [k * I])))[1]

@time identifiability_report = check_identifiability(
    io_equation,
    [a, bi, bw, gam, mu, xi, k]
)

println(identifiability_report)
