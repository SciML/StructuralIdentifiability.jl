# Lee, E. C., Kelly, M. R., Ochocki, B. M., Akinwumi, S. M., Hamre, K. E., Tien, J. H., Eisenberg, M. C.,
# Model distinguishability and inference robustness in mechanisms of cholera transmission and loss of immunity
# https://doi.org/10.1016/j.jtbi.2017.01.032
# Eq. (3)
using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

ode = @ODEmodel(S'(t)=mu - bi * S(t) * I(t) - bw * S(t) * W(t) - mu * S(t) + a * R(t),
                I'(t)=bw * S(t) * W(t) + bi * S(t) * I(t) - (gam + mu) * I(t),
                W'(t)=xi * (I(t) - W(t)),
                R'(t)=gam * I(t) - (mu + a) * R(t),
                y(t)=k * I(t))

@time println(assess_identifiability(ode))
