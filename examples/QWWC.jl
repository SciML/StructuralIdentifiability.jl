include("../io_equation.jl")

#QWWC
logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(
    x'(t) = a * (y(t) - x(t)) + y(t) * z(t),
    y'(t) = b * (x(t) + y(t)) - x(t) * z(t),
    z'(t) = - c * z(t) - d * w(t) + x(t) * y(t),
    w'(t) = e * z(t) - f * w(t) + x(t) * y(t)
)

g = x

@time io_equation = collect(values(find_ioequation(ode, [g])))[1]

println("The number of monomials in the IO-equation is $(length(io_equation))")

@time identifiability_report = check_identifiability(
    io_equation,
    [a, b, c, d, e, f];
    method="GroebnerBasis"
)

println(identifiability_report)
