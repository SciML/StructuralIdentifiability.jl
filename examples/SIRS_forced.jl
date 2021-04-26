# Capistran M., Moreles M., Lara B.
# Parameter Estimation of Some Epidemic Models. The Case of Recurrent Epidemics Caused by Respiratory Syncytial Virus
# doi.org/10.1007/s11538-009-9429-3
# Equations (7)-(11)
using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
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
