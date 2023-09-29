# Globally Identifiable Functions

StructuralIdentifiability.jl provides the function `find_identifiable_functions` that returns the identifiable combinations of the ODE variables.

For example, consider the following model[^1].

```@example
using StructuralIdentifiability # hide
LLW1987 = @ODEmodel(
    x1'(t) = -p1 * x1(t) + p2 * u(t),
    x2'(t) = -p3 * x2(t) + p4 * u(t),
    x3'(t) = -(p1 + p3) * x3(t) + (p4 * x1(t) + p2 * x2(t)) * u(t),
    y1(t) = x3(t)
)
```

Several decades ago, this model was introduced to demonstrate the presence of nontrivial **un**identifiability in nonlinear systems of ODEs. Nowadays, we can automatically find the identifiable combinations of parameters:

```@example
using StructuralIdentifiability # hide
find_identifiable_functions(LLW1987)
```

And even of parameters and states:

```@example
find_identifiable_functions(LLW1987, with_states = true)
```

By default, `find_identifiable_functions` tries to simplify the output functions as much as possible. This feature is useful, for example, in [Model Reparametrization](@ref).

[^1]: > Y. Lecourtier, F. Lamnabhi-Lagarrigue, and E. Walter, [*A method to prove that nonlinear models can be unidentifiable*](https://doi.org/10.1109/CDC.1987.272467), Proceedings of the 26th Conference on Decision and Control, December 1987, 2144-2145;
