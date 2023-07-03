include((@__DIR__) * "/../benchmarks.jl")

ode = Dict(
    :name => "Toy system",
    :ode => @ODEmodel(
        x1'(t) = (a + b)^3 * x1(t) - c * x1(t) * x2(t),
        x2'(t) = -a^3 * b * x2(t) + d^3 * x1(t) * x2(t),
        y1(t) = x1(t)
    ),
    :skip => false,
)
push!(benchmarks, ode)
