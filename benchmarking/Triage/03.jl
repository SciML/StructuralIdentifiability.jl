# https://doi.org/10.3166/ejc.8.343-359
# eq. 8

using StructuralIdentifiability, Nemo

ode = @ODEmodel(
        x1'(t) = -k*x1(t)*x2(t) + k3*x3(t) + k4*x4(t),

        x2'(t) = -k*x1(t)*x2(t) + k3*x3(t) + k4*x4(t),

        x3'(t) = k*x1(t)*x2(t) - (k3 + beta3)*x3(t),

        x4'(t) = beta3*x3(t) - k4*x4(t),

        y1(t) = x1(t)*x2(t)^2,
        y2(t) = x1(t)*x4(t)
    )

# Too hard
# assess_identifiability(ode)

#############################################
# Values from Section 5.1
# Also too hard

param_values = Dict(
    k     => QQ(1, 2),
    k3    => QQ(3),
    k4    => QQ(2),
    beta3 => QQ(1),
)

ode_specialized = set_parameter_values(ode, param_values)
# assess_identifiability(ode_specialized)

#############################################
# Works when reduced with first integrals:
# c1 = x1 - x2,
# c2 = x1 + x3 + x4.

ode_reduced = @ODEmodel(
    x1'(t) = -(1//2)*x1(t)*(x1(t) - c1) + 3*c2 - 3*x1(t) - x4(t),

    x4'(t) = c2 - x1(t) - 3*x4(t),

    y1(t) = x1(t)*(x1(t) - c1)^2,
    y2(t) = x1(t)*x4(t)
)
assess_identifiability(ode_reduced)

#############################################
# Works when calling simplification on y, L^1 y, L^2 y, ...

include(joinpath(@__DIR__, "shared.jl"))

gens, simple_gens = simplify_lie(ode, 4);
println(simple_gens)
