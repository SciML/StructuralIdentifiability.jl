# Creating ODE System

Let us consider the following five-dimensional model with two outputs:

$\begin{cases}x'(t) = lm - d \, x(t) - \beta \, x(t) \, v(t),\\
y'(t) = \beta \, x(t) \, v(t) - a \, y(t),\\
v'(t) = k \, y(t) - u \, v(t),\\
w'(t) = c \, x(t) \, y(t) \, w(t) - c \, q \, y(t) \, w(t) - b \, w(t),\\
z'(t) = c \, q \, y(t) \, w(t) - h \, z(t),\\
y_1(t) = w(t),\\
y_2(t) = z(t)\end{cases}$

This model describes HIV dynamics[^1].

Creating this system in StructuralIdentifiability.jl is easy with the `@ODEmodel` macro:

```@example
using StructuralIdentifiability

ode = @ODEmodel(
    x'(t) = lm - d * x(t) - beta * x(t) * v(t),
    y'(t) = beta * x(t) * v(t) - a * y(t),
    v'(t) = k * y(t) - u * v(t),
    w'(t) = c * x(t) * y(t) * w(t) - c * q * y(t) * w(t) - b * w(t),
    z'(t) = c * q * y(t) * w(t) - h * z(t),
    y1(t) = w(t),
    y2(t) = z(t)
)
```

The resulting `ode` is an instance of the `StructuralIdentifiability.ODE` object.
Most of the functions in StructuralIdentifiability.jl work with `StructuralIdentifiability.ODE`.

For example, we can find the Input-Output equations of this system by passing the `ode` to `find_ioequations`:

```@example
using StructuralIdentifiability # hide
find_ioequations(ode)
```

[^1]: > D. Wodarz, M. Nowak, [*Specific therapy regimes could lead to long-term immunological control of HIV*](https://doi.org/10.1073/pnas.96.25.14464), PNAS December 7, 1999 96 (25) 14464-14469;
