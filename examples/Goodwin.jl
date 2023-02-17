# B.C. Goodwin
# Oscillatory behavior in enzymatic control processes
# https://doi.org/10.1016/0065-2571(65)90067-1

using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

ode = @ODEmodel(x1'(t)=-b * x1(t) + 1 / (c + x4(t)),
                x2'(t)=alpha * x1(t) - beta * x2(t),
                x3'(t)=gama * x2(t) - delta * x3(t),
                x4'(t)=sigma * x4(t) * (gama * x2(t) - delta * x3(t)) / x3(t),
                y(t)=x1(t))

@time println(assess_identifiability(ode))
