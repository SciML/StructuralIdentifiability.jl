include("io_equation.jl")

logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

varnames = vcat(["x$i" for i in 1:7], ["p$i" for i in 1:27], ["u"])

R, Rvars = PolynomialRing(QQ, varnames)
xs = Rvars[1:7]
ps = Rvars[8:(end - 1)]
u = Rvars[end]

ode = ODE(Dict(xs[1] => ps[1] * xs[6] // (ps[3] + xs[6]) - ps[5] * xs[1] // (ps[12] + xs[1]) +
                        ps[26] * xs[7] * u,
               xs[2] => ps[19] * xs[1] - ps[22] * xs[2] + ps[23] * xs[3] -
                        ps[6] * xs[2] // (ps[13] + xs[2]),
               xs[3] => ps[22] * xs[2] - ps[23] * xs[3] - ps[7] * xs[3] // (ps[14] + xs[3]),
               xs[4] => ps[2] * ps[4]^2 // (ps[4]^2 + xs[3]^2) - ps[8] * xs[4] // (ps[15] + xs[4]),
               xs[5] => ps[20] * xs[4] - ps[24] * xs[5] + ps[25] * xs[6] -
                        ps[9] * xs[5] // (ps[16] + xs[5]),
               xs[6] => ps[24] * xs[5] - ps[25] * xs[6] - ps[10] * xs[6] // (ps[17] + xs[6]),
               xs[7] => ps[21] - ps[11] * xs[7] // (ps[18] + xs[7]) - (ps[21] + ps[27] * xs[7]) * u),
          [u])

@time io_equations = find_ioequations(ode, [xs[1], xs[4]])

print(map(length, io_equations), "\n")
