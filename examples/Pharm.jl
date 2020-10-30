include("../io_equation.jl")

logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(
    x0'(t) = a1 * (x1(t) - x0(t)) - (ka * n * x0(t)) / (kc * ka + kc * x2(t) + ka * x0(t)),
    x1'(t) = a2 * (x0(t) - x1(t)),
    x2'(t) = b1 * (x3(t) - x2(t)) - (kc * n * x2(t)) / (kc * ka + kc * x2(t) + ka * x0(t)),
    x3'(t) = b2 * (x2(t) - x3(t))
)

g = x0

@time io_equation = collect(values(find_ioequation(ode, [g])))[1]

println("The number of monomials in the IO-equation is $(length(io_equation))")

@time identifiability_report = check_identifiability(
    io_equation,
    [a1, a2, b1, b2, kc, ka, n];
    method="GroebnerBasis"
)

println(identifiability_report)
