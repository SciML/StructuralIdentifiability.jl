# Lee, E. C., Kelly, M. R., Ochocki, B. M., Akinwumi, S. M., Hamre, K. E., Tien, J. H., Eisenberg, M. C.,
# Model distinguishability and inference robustness in mechanisms of cholera transmission and loss of immunity
# https://doi.org/10.1016/j.jtbi.2017.01.032
# Eq. (3) + extra output
using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(
    x_0'(t) = mu - bi * x_0(t) * x_1(t) - bw * x_0(t) * x_2(t) - mu * x_0(t) + a * x_3(t),
    x_1'(t) = bw * x_0(t) * x_2(t) + bi * x_0(t) * x_1(t) - (gam + mu) * x_1(t),
    x_2'(t) = xi * (x_1(t) - x_2(t)),
    x_3'(t) = gam * x_1(t) - (mu + a) * x_3(t),
    y1(t) = k * x_1(t),
    y2(t) = x_0(t) + x_1(t) + x_3(t)
)

@time println(assess_global_identifiability(ode))
