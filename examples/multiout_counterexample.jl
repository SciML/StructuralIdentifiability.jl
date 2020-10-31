include("../src/io_equation.jl")

logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(
    x1'(t) = (1 + x1(t)^2) // 2,
    x2'(t) = (1 - x1(t)^2) // (1 + x1(t)^2),
    [b]
)

println(ode.parameters)

g = Array{Union{fmpq_mpoly, Generic.Frac{fmpq_mpoly}}, 1}([2 * x1 // (b * (1 + x1^2)), x2])

@time io_equations = find_ioequations(ode, g)

println(io_equations)

@time identifiability_report = check_identifiability(collect(values(io_equations)), ode.parameters)

println(identifiability_report)
