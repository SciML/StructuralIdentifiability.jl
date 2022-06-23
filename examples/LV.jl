# Lotka-Volterra model
using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(x1'(t)=a * x1(t) - b * x1(t) * x2(t) + u(t),
                x2'(t)=-c * x2(t) + d * x1(t) * x2(t),
                y(t)=x1(t))

@time println(assess_identifiability(ode))
