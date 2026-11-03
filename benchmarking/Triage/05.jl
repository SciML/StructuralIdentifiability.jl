# https://doi.org/10.1007/s40313-014-0129-3

using StructuralIdentifiability, Nemo

# eq. 1, two inputs
ode = @ODEmodel(
    x1'(t) = (mu1*S1(t)/(KS1 + S1(t)) - alpha*D)*x1(t),
    x2'(t) = (mu2*S2(t)/(KS2 + S2(t)) - alpha*D)*x2(t),
    S1'(t) = D*(u1(t) - S1(t)) - k1*mu1*S1(t)*x1(t)/(KS1 + S1(t)),
    S2'(t) = D*(u2(t) - S2(t)) + k2*mu1*S1(t)*x1(t)/(KS1 + S1(t)) - k3*mu2*S2(t)*x2(t)/(KS2 + S2(t)),
    y(t) = k6*mu2*S2(t)*x2(t)/(KS2 + S2(t))
)

## too hard
# assess_identifiability(ode)

## works
include(joinpath(@__DIR__, "shared.jl"))

gens, simple_gens = simplify_lie(ode, 5);
println(simple_gens)

# eq. 5, no inputs
ode = @ODEmodel(
    x1'(t) = mu1*S1(t)*x1(t)/(KS1 + S1(t)),
    x2'(t) = mu2*S2(t)*x2(t)/(KS2 + S2(t)),
    S1'(t) = -k1*mu1*S1(t)*x1(t)/(KS1 + S1(t)),
    S2'(t) = k2*mu1*S1(t)*x1(t)/(KS1 + S1(t)) - k3*mu2*S2(t)*x2(t)/(KS2 + S2(t)),
    y(t) = k6*mu2*S2(t)*x2(t)/(KS2 + S2(t))
)

## too hard
# assess_identifiability(ode1)

## doesn't work
include(joinpath(@__DIR__, "shared.jl"))

# gens, simple_gens = simplify_lie(ode, 6);
println(simple_gens)
