# Conradi, C., Shiu, A.,
# Dynamics of post-translational modification systems: recent progress and future directions
# Eq. 3.4
using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

ode = @ODEmodel(x1'(t) = -k1 * x1(t) * x2(t) + k2 * x4(t) + k4 * x6(t),
                x2'(t) = -k1 * x1(t) * x2(t) + k2 * x4(t) + k3 * x4(t),
                x3'(t) = k3 * x4(t) + k5 * x6(t) - k6 * x3(t) * x5(t),
                x4'(t) = k1 * x1(t) * x2(t) - k2 * x4(t) - k3 * x4(t),
                x5'(t) = k4 * x6(t) + k5 * x6(t) - k6 * x3(t) * x5(t),
                x6'(t) = -k4 * x6(t) - k5 * x6(t) + k6 * x3(t) * x5(t),
                y1(t) = x3(t),
                y2(t) = x2(t))

@time println(assess_identifiability(ode))
