include("../src/io_equation.jl")

logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(
    x1'(t) = (1 + x1(t)^2) // 2,
    x2'(t) = (1 - x1(t)^2) // (1 + x1(t)^2),
    y1(t) = 2 * x1 // (b * (1 + x1^2)),
    y2(t) = x2
)

@time io_equations = find_ioequations(ode)

println(io_equations)
