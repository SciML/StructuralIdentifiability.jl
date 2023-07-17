# Demignot, S., D., D., 
# Effect of prosthetic sugar groups on the pharmacokinetics of glucose-oxidase
# https://pubmed.ncbi.nlm.nih.gov/2855567/
using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(
    x0'(t) =
        a1 * (x1(t) - x0(t)) - (ka * n * x0(t)) / (kc * ka + kc * x2(t) + ka * x0(t)),
    x1'(t) = a2 * (x0(t) - x1(t)),
    x2'(t) =
        b1 * (x3(t) - x2(t)) - (kc * n * x2(t)) / (kc * ka + kc * x2(t) + ka * x0(t)),
    x3'(t) = b2 * (x2(t) - x3(t)),
    y1(t) = x0(t)
)

@time println(assess_identifiability(ode))
