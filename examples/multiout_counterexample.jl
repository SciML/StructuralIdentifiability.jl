using Logging

include("../src/StructuralIdentifiability.jl")
using .StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(
    x1'(t) = (1 + x1(t)^2) // 2,
    x2'(t) = (1 - x1(t)^2) // (1 + x1(t)^2),
    y1(t) = 2 * x1 // (b * (1 + x1^2)),
    y2(t) = x2
)

@time println(assess_global_identifiability(ode))
