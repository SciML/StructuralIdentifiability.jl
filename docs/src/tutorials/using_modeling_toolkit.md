# Using `ModelingToolkit.jl` With `StructuralIdentifiability.jl`

In this tutorial, we will cover examples of solving identifiability problems for models defined with the syntax of `ModelingToolkit.jl`.

## Input System

Let us consider the following ODE model with two outputs:

$\begin{cases}
    \dot{S} = -b \, S \, (I + J + q \, A) \, N_{inv},
    \dot{E} = b \, S \, (I + J + q \, A) \, N_{inv} - k \, E,
    \dot{A} = k \, (1 - r) \, E - g_1 \, A,
    \dot{I} = k \, r \, E - (\alpha + g_1) \, I,
    \dot{J} = \alpha \, I - g_2 \, J,
    \dot{C} = \alpha \, I,
    y_1 = C,
    y_2 = N_{inv}
\end{cases}$

The main difference between the input formats in `ModelingToolkit.jl` and `StructuralIdentifiability.jl` is that the output (measured values/functions) must be specified separately in `ModelingToolkit.jl`. In this example, measured quantities are presented by $y_1$, $y_2$.

First, let us define the ODE. We will use `@parameters` and `@variables` macro to define parameters and time-depended functions in the ODE.

```julia
using StructuralIdentifiability, ModelingToolkit

@parameters b q N_inv k r alpha g1 g2
@variables t S(t) E(t) A(t) I(t) J(t) C(t) y1(t) y2(t)
```

The actual ODE will be defined using `ODESystem` structure from `ModelingToolkit.jl`:

```julia
D = Differential(t)

eqs = [
    D(S) ~ -b * S * (I + J + q * A) * N_inv,
    D(E) ~ b * S * (I + J + q * A) * N_inv - k * E,
    D(A) ~ k * (1 - r) * E - g1 * A,
    D(I) ~ k * r * E - (alpha + g1) * I,
    D(J) ~ alpha * I - g2 * J,
    D(C) ~ alpha * I,
]

ode = ODESystem(eqs, t, name = :SEIAJRCmodel)
```

Finally, let us define the array of measured quantities and call the `assess_identifiability` function. This is the main function that determines local/global identifiability properties of each parameter and state. We will use the probability of correctness $p=0.99$.

For `ModelingToolkit.jl`, both `assess_identifiability` and `assess_local_identifiability` functions accept keyword arguments: 

* `measured_quantities`, also called "output functions" in identifiability literature; these are crucial for answering identifiability questions.
* `p`, probability of correctness. This value equals 0.99 by default.
* `funcs_to_check`, functions of parameters of which we wish to check identifiability.

```julia
measured_quantities = [y1 ~ C, y2 ~ N_inv]
@time global_id = assess_identifiability(ode, measured_quantities=measured_quantities)
```

Let us put all of the code above together:

```@repl
using StructuralIdentifiability, ModelingToolkit

@parameters b q N_inv k r alpha g1 g2
@variables t S(t) E(t) A(t) I(t) J(t) C(t) y1(t) y2(t)

D = Differential(t)

eqs = [
    D(S) ~ -b * S * (I + J + q * A) * N_inv,
    D(E) ~ b * S * (I + J + q * A) * N_inv - k * E,
    D(A) ~ k * (1 - r) * E - g1 * A,
    D(I) ~ k * r * E - (alpha + g1) * I,
    D(J) ~ alpha * I - g2 * J,
    D(C) ~ alpha * I,
]

ode = ODESystem(eqs, t, name = :SEIAJRCmodel)

measured_quantities = [y1 ~ C, y2 ~ N_inv]
@time global_id = assess_identifiability(ode, measured_quantities=measured_quantities)
```
<!-- Dict{Num, Symbol} with 8 entries:
  k     => :globally
  b     => :globally
  alpha => :globally
  g1    => :globally
  g2    => :globally
  r     => :nonidentifiable
  q     => :nonidentifiable
  N_inv => :globally -->
  
<!-- Indeed, notice how much quicker we obtained the result with 99% correctness guarantee! This illustrates the fact that you may sometimes sacrifice probability slightly to get results much faster. -->

[^1]:
    > K. Roosa and G. Chowell. [*Assessing parameter identifiability in compartmental dynamic models using a computational approach: application to infectious disease transmission models*](https://doi.org/10.1186/s12976-018-0097-6), Theor Biol Med Model 16, 1 (2019)