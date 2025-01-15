# Reparametrizations

## Overview

Once one has found that not all parameters and/or states of the model at hand are identifiable, one natural desire is to
reparametrize the model into a one with better identifiability properties.
`StructuralIdentifiability` offers such a functionality via the function `reparametrize_global`.
It takes as input an ODE model and produces its transformation into another model with the same
input-output behaviour but with the states and parameters being globally identifiable.
Note that, in general, such a transformation may not exist in the class of rational models,
so sometimes the function returns an ODE not on the whole affine space but on a manifold.

More precisely, the function returns a dictionary with three keys:

  - `:new_vars` is a dictionary which maps the new parameters and new states into the formulas expressing them in terms of the original parameters and states;

  - `:new_ode` is the ODE satisfied by these new states (and the expression of the output in terms of the new states);
  - `:implicit_relations` is a list of algebraic relations between the new states and parameters. Being nonempty, this is exactly the list of equations defining the manifold, on which the new ODE model is defined. In many interesting  cases, however, this list is empty meaning that the new ODE is a standard rational ODE model.

## Example: SEUIR model

Consider a SEUIR epidemiological model from[^1]:

$\begin{cases}
S(t)' = -\beta \frac{(U(t) + I(t))S(t)}{N},\\
E(t)' = \beta \frac{(U(t) + I(t))S(t)}{N} - \gamma E(t),\\
U(t)' = (1 - \alpha) \gamma E(t) - \delta U(t),\\
I(t)' = \alpha \gamma E(t) - \delta I(t),\\
R(t)' = \delta(U(t) + I(t)),\\
y(t) = I(t)
\end{cases}$

In this model `S` is, as usually, the number of susceptible people, `E` is the number of people exposed to virus but not yet infected
(as in a simple SEIR model[^1]), and `I` and `U` correspond to number of infected people who report the infection and who do not, respectively.
We define the model but omit `R` compartment since it does not affect the output dynamics:

```@example seuir
using StructuralIdentifiability

ode = @ODEmodel(
    S'(t) = -b * (U(t) + I(t)) * S(t) / N,
    E'(t) = b * (U(t) + I(t)) * S(t) / N - g * E(t),
    U'(t) = (1 - a) * g * E(t) - d * U(t),
    I'(t) = a * g * E(t) - d * I(t),
    y(t) = I(t)
)
```

Majority of the states and parameters are not identifiable in this case:

```@example seuir
assess_identifiability(ode)
```

Let us attempt to reparametrize the model, and print new variables:

```@example seuir
reparam = reparametrize_global(ode)
@assert isempty(reparam[:implicit_relations]) # checking that the result is an ODE on the whole space, not on a manifold
reparam[:new_vars]
```

In these new variables and parameters, the original ODE can be rewritten as follows:

```@example seuir
reparam[:new_ode]
```

In order to analyze this result, let us give more interpretable names to the new variables and parameters:

$I := I, \; \widetilde{E} := \alpha E, \widetilde{S} := \alpha, \; \widetilde{I} := \alpha (I + U), \; \gamma := \gamma,\;\delta := \delta,\;\widetilde{N} := \frac{\alpha N}{\beta}$

Then the reparametrize system becomes

$\begin{cases}
\widetilde{S}'(t) = -\widetilde{S}(t) \widetilde{I}(t) / \widetilde{N},\\
\widetilde{E}'(t) = \widetilde{S}(t) \widetilde{I}(t) / \widetilde{N} - \gamma \widetilde{E}(t),\\
\widetilde{I}'(t) = -\delta \widetilde{I}(t) + \gamma\widetilde{E}(t),\\
I'(t) = \gamma\widetilde{E}(t) - \delta I(t),\\
y(t) = I(t)
\end{cases}$

This reparametrization not only reduces the dimension of the parameter space from 5 to 3 but reveals interesting structural properties of the model:

  - The first three equations form a self-contained model which is equivalent to a simple SEIR model, so the model gets "decomposed";

  - New variables $\widetilde{S}$, $\widetilde{E}$, $\widetilde{I}$ are obtained from $S$, $E$, and $I$ by scaling by $\alpha$ which is the ratio of people who report being infected. One can interpret this as there is a part of population who would report infection and the other part who would not. Ultimately, we can model only the ones who would as this is mainly they who contribute to the output.

Finally, we can check that the new model is indeed globally identifiable:

```@example seuir
assess_identifiability(reparam[:new_ode])
```

[^1]: > T. Sauer, T. Berry, D. Ebeigbe, M. Norton, A. Whalen, S. Schiff, [*Identifiability of infection model parameters early in an epidemic*](https://doi.org/10.1137/20m1353289), SIAM Journal on Control and Optimization, 2022;
