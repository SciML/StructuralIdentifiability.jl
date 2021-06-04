using StructuralIdentifiability
using StructuralIdentifiability: numerical_identifiability

using HomotopyContinuation

ode = @ODEmodel(
    x'(t) = x(t) * a,
    y(t) = x(t) + b^2
)

#ode = @ODEmodel(
#    S'(t) = mu - bi * S(t) * I(t) - bw * S(t) * W(t) - mu * S(t) + a * R(t),
#    I'(t) = bw * S(t) * W(t) + bi * S(t) * I(t) - (gam + mu) * I(t),
#    W'(t) = xi * (I(t) - W(t)),
#    R'(t) = gam * I(t) - (mu + a) * R(t),
#    y(t) = k * I(t)
#)
#ode = @ODEmodel(
#    S'(t) = 1 - 1//3 * S(t)^2 * I(t) - 2 * S(t) * W(t) - 1 * S(t) + 4 * R(t),
#    I'(t) = 2 * S(t) * W(t) + 1//3 * S(t) * I(t) - (-2 + 1) * I(t),
#    W'(t) = 1 * (I(t)^4 - W(t)),
#    R'(t) = -2 * I(t) - (1 + 4) * R(t),
#    #bi'(t) = 0,
#    #k'(t) = 0,
#    y(t) = 1//7 * I(t)
#)

println(numerical_identifiability(ode))
