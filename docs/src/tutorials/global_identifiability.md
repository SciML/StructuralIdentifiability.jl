# Global Identifiability of Differential Models

In this tutorial, let us cover an example problem of querying the ODE for globally identifiable parameters.

## Input System

Let us consider the following four-dimensional model with two outputs:

$\begin{cases}x'(t) = lm - d \, x(t) - \beta \, x(t) \, v(t),\\
    y'(t) = \beta \, x(t) \, v(t) - a \, y(t),\\
    v'(t) = k \, y(t) - u \, v(t),\\
    w'(t) = c \, x(t) \, y(t) \, w(t) - c \, q \, y(t) \, w(t) - b \, w(t),\\
    z'(t) = c \, q \, y(t) \, w(t) - h \, z(t),\\
    y_1(t) = w(t),\\
    y_2(t) = z(t)\end{cases}$

This model describes HIV dynamics[^1]. Let us run a global identifiability check on this model to get the result with probability of correctness being `p=0.99`. To do this, we will use `assess_identifiability(ode, p)` function.

Global identifiability needs information about local identifiability first, hence the function we chose here will take care of that extra step for us.

```@repl
using BenchmarkTools
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
@btime global_id = assess_identifiability(ode, 0.99)
```

We also note that it's usually inexpensive to obtain the result with higher probability of correctness. Take, for example, the previous system with `p=0.9999`

```@repl
using BenchmarkTools
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
@btime global_id = assess_identifiability(ode, 0.9999)
```

[^1]:
    > D. Wodarz, M. Nowak, [*Specific therapy regimes could lead to long-term immunological control of HIV*](https://doi.org/10.1073/pnas.96.25.14464), PNAS December 7, 1999 96 (25) 14464-14469;
