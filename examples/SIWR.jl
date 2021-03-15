using Logging

include("../src/StructuralIdentifiability.jl")
using .StructuralIdentifiability

# SIWR Cholera model
logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(
    S'(t) = mu - bi * S(t) * I(t) - bw * S(t) * W(t) - mu * S(t) + a * R(t),
    I'(t) = bw * S(t) * W(t) + bi * S(t) * I(t) - (gam + mu) * I(t),
    W'(t) = xi * (I(t) - W(t)),
    R'(t) = gam * I(t) - (mu + a) * R(t),
    y(t) = k * I(t)
)

@time println(assess_global_identifiability(ode))
