using Logging

include("../src/StructuralIdentifiability.jl")
using .StructuralIdentifiability

# SIWR Cholera model
logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(
  s'(t) = mu - mu * s(t) - b0 * (1 + b1 * x1(t)) * i(t) * s(t) + g * r(t),
  i'(t) = b0 * (1 + b1 * x1(t)) * i(t) * s(t) - (nu + mu) * i(t),
  r'(t) = nu * i(t) - (mu + g) * r(t),
  x1'(t) = -M * x2(t),
  x2'(t) = M * x1(t),
  y1(t) = i(t),
  y2(t) = r(t)
)

@time println(assess_identifiability(ode))
