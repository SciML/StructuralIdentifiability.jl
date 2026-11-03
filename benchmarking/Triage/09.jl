# https://doi.org/10.1371/journal.pone.0027755
# eq. 25

using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = n1*x6(t)/(g1+x6(t)) - m1*x1(t)/(k1+x1(t)) + q1*x7(t)*u(t),
    x2'(t) = p1*x1(t) - r1*x2(t) + r2*x3(t) - m2*x2(t)/(k2+x2(t)),
    x3'(t) = r1*x2(t) - r2*x3(t) - m3*x3(t)/(k3+x3(t)),
    x4'(t) = n2*g2^2/(g2^2+x3(t)^2) - m4*x4(t)/(k4+x4(t)),
    x5'(t) = p2*x4(t) - r3*x5(t) + r4*x6(t) - m5*x5(t)/(k5+x5(t)),
    x6'(t) = r3*x5(t) - r4*x6(t) - m6*x6(t)/(k6+x6(t)),
    x7'(t) = p3 - m7*x7(t)/(k7+x7(t)) - (p3+q2*x7(t))*u(t),
    y1(t) = x1(t),
    y2(t) = x4(t)
)

### Seems hard
# assess_identifiability(ode)

### Seems hard
include(joinpath(@__DIR__, "shared.jl"))
# gens, simple_gens = simplify_lie(ode, 4);
# println(simple_gens)

# eq. 23
ode = @ODEmodel(
    x1'(t) = -k1*x1(t)*u1(t)/(x1(t)+kM),
    x2'(t) = k1*x1(t)*u1(t)/(x1(t)+kM) - k2*x2(t)*u2(t)/(x2(t)+kM),
    x3'(t) = k2*x2(t)*u2(t)/(x2(t)+kM) - k3*x3(t)*u3(t)/(x3(t)+kM),
    x4'(t) = k2*x2(t)*u2(t)/(x2(t)+kM) + k3*x3(t)*u3(t)/(x3(t)+kM) - k4*x4(t)*u4(t)/(x4(t)+kM),
    x5'(t) = k4*x4(t)*u4(t)/(x4(t)+kM),
    y1(t) = x1(t),
    y2(t) = x2(t),
    y3(t) = x3(t),
    y4(t) = x4(t),
    y5(t) = x5(t)
)

# fully identifiable
assess_identifiability(ode)
