include("../src/io_equation.jl")

logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(
    x1'(t) = k_prod  - k_deg * x1(t) - k1 * x1(t) * u(t),
    x2'(t) = -k3 * x2(t) - k_deg * x2(t) a2 * x2(t) * x10(t) + t1 * x4(t) - a3 * x2(t) * x13(t) + t2 * x5(t) + (k1 * x1(t) - k2 * x2(t) * x8(t)) * u(t),
    x3'(t) = k3 * x2(t) - k_deg * x3(t) + k2 * x2(t) * x8(t) * u(t),
    x4'(t) = a2 * x2(t) * x10(t) - t1 * x4(t),
    x5'(t) = a3 *  x2(t) * x13(t) - t2 *  x5(t),
    x6'(t) = c6a * x13(t) - a1 * x6(t) * x10(t) - t2 * x5(t) - i1 * x6(t),
    x7'(t) = i1 * kv * x6(t) - a1 * x11(t) *  x7(t),
    x8'(t) = c4 * x9(t) - c5 * x8(t),
    x9'(t) = c2 - c1 * x7(t) - c3 * x9(t),
    x10'(t) = -a2 * x2(t) * x10(t) - a1 * x10(t) * x6(t) + c4a * x12(t) - c5a * x10(t) - i1a * x10(t) + e1a * x11(t),
    x11'(t) = -a1 * x11(t) * x7(t) + i1a * kv * x10(t) - e1a * kv * x11(t),
    x12'(t) = c2a + c1a * x7(t) - c3a * x12(t),
    x13'(t) = a1 * x10(t) * x6(t) - c6a * x13(t) - a3 * x2(t) * x13(t) + e2a * x14(t),
    x14'(t) = a1 * x11(t) * x7(t) - e2a * kv * x14(t),
    x15'(t) = c2c + c1c * x7(t) - c3c * x15(t)
)

@time io_equations = find_ioequations(ode, Array{Union{fmpq_mpoly, Generic.Frac{fmpq_mpoly}}, 1}([x7, x10 + x13, x9, x1 + x2 + x3, x2, x12]))

eq_list = collect(values(io_equations))

println(map(length, eq_list), "\n")

@time identifiability_report = check_identifiability(eq_list, ode.parameters)

println(identifiability_report)
