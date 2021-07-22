# Local Identifiability of Differential Models

In this tutorial, let us cover an example problem of querying the ODE for globally identifiable parameters.

## Input System

Let us consider the following model with two outputs and a system of four ordinary differential equations:

$\begin{cases}x'(t) = lm - d * x(t) - beta * x(t) * v(t),\\
    y'(t) = beta * x(t) * v(t) - a * y(t),\\
    v'(t) = k * y(t) - u * v(t),\\
    w'(t) = c * x(t) * y(t) * w(t) - c * q * y(t) * w(t) - b * w(t),\\
    z'(t) = c * q * y(t) * w(t) - h * z(t),\\
    y_1(t) = w(t),\\
    y_2(t) = z(t)\end{cases}$

This model comes from [[1]](#hiv), describing HIV dynamics. Let us run a global identifiability check on this model to get the result with probability of correctness being `p=0.999`.

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
@time global_id = assess_global_identifiability(ode, 0.999)
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
@time global_id = assess_global_identifiability(ode, 0.99)
```

Indeed, notice how much quicker we obtained the result with 99% correctness guarantee! This illustrates the fact that you may sometimes sacrifice probability slightly to get results much faster.

## References

<a id="hiv">[1]</a> Wodarz, D., Nowak, M.,
[*Specific therapy regimes could lead to long-term immunological control of HIV*](https://doi.org/10.1073/pnas.96.25.14464),
Proceedings of the National Academy of Sciences Dec 1999, 96 (25)

