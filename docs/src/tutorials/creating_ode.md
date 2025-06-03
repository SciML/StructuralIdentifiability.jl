# Creating ODE System

Most of the algorithms in `StructuralIdentifiability.jl` take as input system of ordinary differential equations (ODEs)
in the state space form, that is:

$\begin{cases}
\mathbf{x}'(t) = \mathbf{f}(\mathbf{x}(t), \mathbf{p}, \mathbf{u}(t)),\\
\mathbf{y}(t) = \mathbf{g}(\mathbf{x}(t), \mathbf{p}, \mathbf{u(t)}),
\end{cases}$

which involves

  - a vector $\mathbf{x}(t)$ of the state variables of the system,

  - a vector $\mathbf{u}(t)$ of external inputs,
  - a vector $\mathbf{p}$ of scalar parameters,
  - a vector $\mathbf{y}(t)$ of outputs (i.e., observations),
  - and vectors of rational functions $\mathbf{f}$ and $\mathbf{g}$ (for discussion of the non-rational case, see this [issue](https://github.com/SciML/StructuralIdentifiability.jl/issues/144)).

In the standard setup, inputs and outputs are assumed to be known, and the goal is to assess
**identifiability** of parameters and/or states from the input-output data.
In the case of states, this property is also called **observability**.

There are two ways to define such a system to be processed using `StructuralIdentifiability.jl`.
We will demonstrate them using the following example system
(Wright's population model of two mutualist species with control[^1]):

$\begin{cases}
x_1'(t) = r_1 x_1(t)(1 - c_1 x_1(t)) + \frac{\beta_1 x_1(t)x_2(t)}{\chi_1 + x_2(t)} + u(t),\\
x_2'(t) = r_2 x_2(t)(1 - c_2 x_2(t)) + \frac{\beta_2 x_1(t)x_2(t)}{\chi_2 + x_1(t)},\\
y(t) = x_1(t).
\end{cases}$

## Defining the model using `@ODEmodel` macro

One way to define the model is to use the `@ODEmodel` macro provided by the `StructuralIdentifiability.jl` package.

```@example 1
using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) =
        r1 * x1(t) * (1 - c1 * x1(t)) + beta1 * x1(t) * x2(t) / (chi1 + x2(t)) + u(t),
    x2'(t) = r2 * x2(t) * (1 - c2 * x2(t)) + beta2 * x1(t) * x2(t) / (chi2 + x1(t)),
    y(t) = x1(t)
)
```

Then one can, for example, assess identifiability of the parameters and states by

```@example 1
assess_identifiability(ode)
```

## Defining using `ModelingToolkit`

`StructuralIdentifiability` has an extension `ModelingToolkitSIExt` which allows to use `System` from `ModelingToolkit` to describe
a model. The extension is loaded automatically once `ModelingToolkit` is loaded via `using ModelingToolkit`.
In this case, one should encode the equations for the states as `System` and specify the outputs separately.
In order to do this, we first introduce all functions and scalars:

```@example 2; continued = true
using StructuralIdentifiability, ModelingToolkit

@parameters r1, r2, c1, c2, beta1, beta2, chi1, chi2
@variables t, x1(t), x2(t), y(t), u(t)

D = Differential(t)
```

And then defined the system and separately the outputs:

```@example 2
eqs = [
    D(x1) ~ r1 * x1 * (1 - c1 * x1) + beta1 * x1 * x2 / (chi1 + x2) + u,
    D(x2) ~ r2 * x2 * (1 - c2 * x2) + beta2 * x1 * x2 / (chi2 + x1),
]

measured_quantities = [y ~ x1]

ode_mtk = System(eqs, t, name = :mutualist)
```

Then, for example, the identifiability of parameters and states can be assessed as follows:

```@example 2
assess_identifiability(ode_mtk, measured_quantities = measured_quantities)
```

[^1]: > D. H. Wright, [*A Simple, Stable Model of Mutualism Incorporating Handling Time*](https://doi.org/10.1086/285003), The American Naturalist, 1989, 134(4).
