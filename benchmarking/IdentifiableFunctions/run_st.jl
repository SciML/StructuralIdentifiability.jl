#! format: off
using StructuralIdentifiability
using Logging
Base.global_logger(ConsoleLogger(Logging.Info))

St = StructuralIdentifiability.@ODEmodel(
    S'(t) = r * S(t) - (e + a * W(t)) * S(t) - d * W(t) * S(t) + g * R(t),
    R'(t) = rR * R(t) + (e + a * W(t)) * S(t) - dr * W(t) * R(t) - g * R(t),
    W'(t) = Dd * (T - W(t)),
    y1(t) = S(t) + R(t),
    y2(t) = T
)

hiv = StructuralIdentifiability.@ODEmodel(
    w'(t) = -b * w(t) + c * w(t) * x(t) * y(t) - c * w(t) * q * y(t),
    v'(t) = k * y(t) - v(t) * u,
    x'(t) = lm - x(t) * d - x(t) * v(t) * beta,
    z'(t) = c * w(t) * q * y(t) - h * z(t),
    y'(t) = x(t) * v(t) * beta - a * y(t),
    y2(t) = z(t),
    y1(t) = w(t)
)

funcs0 = StructuralIdentifiability.find_identifiable_functions(
    hiv, 
    strategy = (:hybrid,)
)
