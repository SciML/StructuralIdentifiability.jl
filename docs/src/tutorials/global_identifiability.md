# Global Identifiability of Differential Models

In this tutorial, let us cover an example problem of querying the ODE for globally identifiable parameters.

## Input System

Let us consider the following four-dimensional model with two outputs:

$\begin{cases}x'(t) = lm - d * x(t) - beta * x(t) * v(t),\\
    y'(t) = beta * x(t) * v(t) - a * y(t),\\
    v'(t) = k * y(t) - u * v(t),\\
    w'(t) = c * x(t) * y(t) * w(t) - c * q * y(t) * w(t) - b * w(t),\\
    z'(t) = c * q * y(t) * w(t) - h * z(t),\\
    y_1(t) = w(t),\\
    y_2(t) = z(t)\end{cases}$

This model describes HIV dynamics[^1]. Let us run a global identifiability check on this model to get the result with probability of correctness being `p=0.999`. To do this, we will use `assess_identifiability(ode, p)` function.

Global identifiability needs information about local identifiability first, hence the function we chose here will take care of that extra step for us.

```@repl
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
@time global_id = assess_identifiability(ode, 0.999)
```

Now let us compare the same system but with probability being `p=0.99`. We will see a reduction in runtime:

```@repl
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
@time global_id = assess_identifiability(ode, 0.99)
```

Indeed, notice how much quicker we obtained the result with 99% correctness guarantee! This illustrates the fact that you may sometimes sacrifice probability slightly to get results much faster.

[^1]:
    > D. Wodarz, M. Nowak, [*Specific therapy regimes could lead to long-term immunological control of HIV*](https://doi.org/10.1073/pnas.96.25.14464), PNAS December 7, 1999 96 (25) 14464-14469;