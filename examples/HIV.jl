# Wodarz, D., Nowak, M.
# Specific therapy regimes could lead to long-term immunological control of HIV
# https://doi.org/10.1073/pnas.96.25.14464
# Page 1
using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

ode = @ODEmodel(
    x'(t) = lm - d * x(t) - beta * x(t) * v(t),
    y'(t) = beta * x(t) * v(t) - a * y(t),
    v'(t) = k * y(t) - u * v(t),
    w'(t) = c * x(t) * y(t) * w(t) - c * q * y(t) * w(t) - b * w(t),
    z'(t) = c * q * y(t) * w(t) - h * z(t),
    y1(t) = w(t),
    y2(t) = z(t)
)

@time println(assess_identifiability(ode))
