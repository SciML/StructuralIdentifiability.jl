#=
1. We consider the system
    x1' = -M * x2,
    x2' = M * x1.

The observed output is
    y = x1.

The state x2 is not observable.

Parameter M itself is not identifiable.
However, the square of M is identifiable.
Also, the product of M and x2 is identifiable. 
=#
import Pkg;
Pkg.activate(@__DIR__);

using DiffEqParamEstim, DifferentialEquations, RecursiveArrayTools, Plots
using Optimization, ForwardDiff, OptimizationOptimJL, OptimizationBBO
using StructuralIdentifiability

# 2. Define the system
function f(du, u, p, t)
    du[1] = -p[1] * u[2]
    du[2] = p[1] * u[1]
end

# 3. Generate the time-series x1(t), x2(t) that solve the above system for M = 1.5
p = [1.5]
u0 = [1.0; 1.0]
tspan = (0.0, 10.0)
prob = ODEProblem(f, u0, tspan, p)
sol = solve(prob, Tsit5())
t = collect(range(0, stop = 10, length = 200))
data_x1 = VectorOfArray([sol(t[i]) for i in 1:length(t)])[1, :]

# 4. Estimate M from x1(t) and some initial approximation. 
function build_cost_function(
    prob,
    alg,
    loss,
    t;
    prob_generator = DiffEqParamEstim.STANDARD_PROB_GENERATOR,
)
    # adapted from build_loss_objective in DiffEqParamEstim.jl
    cost_function = function (p, _ = nothing)
        tmp_prob = prob_generator(prob, p)
        sol = solve(
            tmp_prob,
            alg;
            saveat = t,
            save_everystep = false,
            dense = false,
            verbose = false,
            maxiters = 10_000,
        )
        failure = !SciMLBase.successful_retcode(sol.retcode)
        failure && return Inf
        sol_x1 = VectorOfArray([sol(t[i]) for i in 1:length(t)])[1, :]
        loss(sol_x1)
    end

    OptimizationFunction(cost_function, Optimization.AutoForwardDiff())
end
# Take M = 1.4 as the initial approximation.
p0 = [1.4]
l2_loss = estimated -> sum((data_x1 - estimated) .^ 2)
cost_function = build_cost_function(prob, Tsit5(), l2_loss, t)
# Sanity check
@assert cost_function(p[1]) < 1e-2 && cost_function(p0[1]) > 1.0

optprob = Optimization.OptimizationProblem(cost_function, p0)
optsol = solve(optprob, BFGS())
@show optsol
@assert abs(p[1] - optsol[1]) < 1e-2

# 5. Plot the landscape of the loss function
vals = -2.0:0.01:2.0
plot(
    vals,
    [cost_function(i) for i in vals],
    yscale = :log10,
    xaxis = "Parameter \$M\$",
    yaxis = "Cost",
    title = "L2 cost function",
    lw = 3,
)
savefig((@__FILE__)[1:(end - 3)] * "-l2.pdf")
