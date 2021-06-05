using Printf
using StructuralIdentifiability
using StructuralIdentifiability: numerical_identifiability

using HomotopyContinuation

# Since the total number of solutions is small, in some runs, for this model
# one gets the result "identifiable"
#ode = @ODEmodel(
#    x'(t) = x(t) * a,
#    y(t) = x(t) + b^2
#)

ode = @ODEmodel(
    x1'(t) = -(a01 + a21) * x1(t) + a12 * x2(t) + u(t),
    x2'(t) = a21 * x1(t) - a12 * x2(t),
    y(t) = x2(t)
)

# Takes a couple of hours!
#ode = @ODEmodel(
#    S'(t) = mu - bi * S(t) * I(t) - bw * S(t) * W(t) - mu * S(t) + a * R(t),
#    I'(t) = bw * S(t) * W(t) + bi * S(t) * I(t) - (gam + mu) * I(t),
#    W'(t) = xi * (I(t) - W(t)),
#    R'(t) = gam * I(t) - (mu + a) * R(t),
#    y(t) = k * I(t)
#)

println(numerical_identifiability(ode))
