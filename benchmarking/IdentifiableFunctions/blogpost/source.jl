import Pkg;
Pkg.activate(@__DIR__);

using DiffEqParamEstim, DifferentialEquations, RecursiveArrayTools, Plots
using Optimization, ForwardDiff, OptimizationOptimJL, OptimizationBBO
using StructuralIdentifiability

#=
d/dt mice (t) = p1 * mice - p2 * mice * owl
d/dt owl  (t) = -p3 * owl + p4 * mice * owl
=#

######################
######################

#!format: off
ode = @ODEmodel(
    🐁'(t) = p1^2 * 🐁 - 🐁 * 🦉, 
    🦉'(t) = -p3 * 🦉 + p4 * 🐁 * 🦉,
    y(t) = 🐁
)
ode = @ODEmodel(
    x1'(t) = -M * x2(t),
	x2'(t) = M * x1(t),
    y(t) = x1
)

find_identifiable_functions(ode, with_states = true)

ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = (a + b)*x1(t) + a*b,
    y1(t) = x1(t),
)

find_identifiable_functions(ode)
#!format: on

######################
######################

function f(du, u, p, t)
    # du[1] = dx = p[1] * u[1] - u[1] * u[2]
    # du[2] = dy = -3 * u[2] + u[1] * u[2]
    du[1] = dx = p[1]^2 * u[1] - u[1] * u[2]
    du[2] = dy = -u[2] + u[1] * u[2]
end

u0 = [1.0; 1.0]
tspan = (0.0, 10.0)
p = [1.5]
prob = ODEProblem(f, u0, tspan, p)

######################
######################

sol = solve(prob, Tsit5())
t = collect(range(0, stop = 10, length = 200))
using RecursiveArrayTools # for VectorOfArray
randomized = VectorOfArray([(sol(t[i]) + 0.01randn(2)) for i in 1:length(t)])
data = convert(Array, randomized)

######################
######################

newprob = remake(prob, p = [1.42])
newsol = solve(newprob, Tsit5())
plot(sol)
plot!(newsol)

######################
######################

@edit L2Loss(t, data)

cost_function = build_loss_objective(
    prob,
    Tsit5(),
    L2Loss(t, data),
    Optimization.AutoForwardDiff(),
    maxiters = 10000,
    verbose = false,
)

######################
######################

vals = -2.0:0.1:2.0
plot(
    vals,
    [cost_function(i) for i in vals],
    yscale = :log10,
    xaxis = "Parameter",
    yaxis = "Cost",
    title = "1-Parameter Cost Function",
    lw = 3,
)

######################
######################

optprob = Optimization.OptimizationProblem(cost_function, [1.4])
optsol = solve(optprob, BFGS())

######################
######################

newprob = remake(prob, p = optsol.u)
newsol = solve(newprob, Tsit5())
plot(sol)
plot!(newsol)
