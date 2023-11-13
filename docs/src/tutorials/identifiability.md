# Identifiability of Differential Models (Local and Global)

Recall that we consider ODE models in the state-space form

$\begin{cases}
\mathbf{x}'(t) = \mathbf{f}(\mathbf{x}(t), \mathbf{p}, \mathbf{u}(t)),\\
\mathbf{y}(t) = \mathbf{g}(\mathbf{x}(t), \mathbf{p}, \mathbf{u(t)}),
\end{cases}$

where $\mathbf{x}(t), \mathbf{y}(t)$, and $\mathbf{u}(t)$ are time-dependent states, outputs, and inputs, respectively,
and $\mathbf{p}$ are scalar parameters.
We will call that a parameter or a states (or a function of them) is **identifiable** if its value can be recovered from
time series for inputs and outputs.
Typically, two types of identifiability are distinguished

  - **local** identifiability: the value can be recovered up to finitely many options;

  - **global** identifiability: the value can be recovered uniquely.

Note that in the case of states, term **observability** it typically used. We will use **identifiability** for both
states and parameters for brevity and uniformity.
While the notion of global identifiability is much more precise, assessing local identifiability is typically much faster,
and can be performed for the models whose global identifiability analysis is out of reach.

## Local identifiability

We consider the Lotka-Volterra model:

$\begin{cases}
x_1'(t) = a x_1(t) - b x_1(t) x_2(t) + u(t),\\
x_2'(t) = -c x_2(t) + d x_1(t) x_2(t),\\
y(t) = x_1(t)
\end{cases}$

The local identifiability of all parameters and states in this model can be assessed as follows

```@example local
using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = a * x1(t) - b * x1(t) * x2(t) + u(t),
    x2'(t) = -c * x2(t) + d * x1(t) * x2(t),
    y(t) = x1(t)
)

assess_local_identifiability(ode)
```

We see that $x_1(t)$ is locally identifiable (no surprises, this is an output), $a, c,$ and $d$ are identifiable as well.
On the other hand, $x_2(t)$ and $b$ are nonidentifiable. This can be explained by the following scaling symmetry

$x_2(t) \to \lambda x_2(t), \quad b \to \frac{b}{\lambda}$

which preserves input and output for every nonzero $\lambda$.
The algorithm behind this call is the one due to Sedoglavic[^1].

Function `assess_local_identifiability` has several optional parameters

  - `funcs_to_check` a list of specific functions of parameters and states to check identifiability for (see an example below).
    If not provided, the identifiability is assessed for all parameters and states.

  - `p` (default $0.99$) is the probability of correctness. The algorithm can, in theory, produce wrong result, but the probability that it is correct
    is guaranteed to be at least `p`. However, the probability bounds we use are quite conservative, so the actual probability of correctness is
    likely to be much higher.

 - `type` (default `:SE`). By default, the algorithm checks the standard single-experiment identifiability. If one sets `type = :ME`, then the algorithm
    checks multi-experiment identifiability, that is, identifiability from several experiments with independent initial conditions (the algorithm from [^2] is used).

 - `loglevel` (default `Logging.Info`). The minimal level of logging messages to be displayed. Available options: `Logging.Debug`, 
   `Logging.Info`, `Logging.Warn`, and `Logging.Error`.

Note that the scaling symmetry given above suggests that $b x_2(t)$ may in fact be identifiable. This can be checked using `funcs_to_check` parameter:

```@example local
assess_local_identifiability(ode, funcs_to_check = [b * x2])
```

Indeed!

## Global identifiability

One can obtain more refined information about a model using `assess_identifiability` function.
We will showcase it using the Goodwin oscillator model [^3].

```@example global
using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = -b * x1(t) + 1 / (c + x4(t)),
    x2'(t) = alpha * x1(t) - beta * x2(t),
    x3'(t) = gama * x2(t) - delta * x3(t),
    x4'(t) = sigma * x4(t) * (gama * x2(t) - delta * x3(t)) / x3(t),
    y(t) = x1(t)
)

assess_identifiability(ode)
```

As a result, each parameter/state is assigned one of the labels `:globally` (globally identifiable), `:locally` (locally but not globally identifiable),
or `:nonidentifiable` (not identifiable, even locally).
The algorithm behind this computation follows [^4].

Similarly to `assess_local_identifiability`, this function has optional parameters:

  - `funcs_to_check` a list of specific functions of parameters and states to check identifiability for (see an example below).
    If not provided, the identifiability is assessed for all parameters and states. Note that the computations for states may be
    more involved than for the parameters, so one may want to call the function with `funcs_to_check = ode.parameters` if the
    call `assess_identifiability(ode)` takes too long.

  - `p` (default $0.99$) is the probability of correctness. Same story as above: the probability estimates are very conservative, so the actual
    error probability is much lower than 1%.
    Also, currently, the probability of correctness does not include the probability of correctness of the modular reconstruction for Groebner bases.
    This probability is ensured by an additional check modulo a large prime, and can be neglected for practical purposes.

  - `loglevel` (default `Logging.Info`). The minimal level of logging messages to be displayed. Available options: `Logging.Debug`, 
    `Logging.Info`, `Logging.Warn`, and `Logging.Error`.

Using `funcs_to_check` parameter, one can further inverstigate the nature of the lack of identifiability in the model at hand.
For example, for the Goodwin oscillator, we can check if `beta + delta` and `beta * delta` are identifiable:

```@example global
assess_identifiability(ode, funcs_to_check = [beta + delta, beta * delta])
```

And we see that they indeed are. This means, in particular, that the reason why `beta` and `delta` are not identifiable is because their values
can be exchanged. One may wonder how could we guess these functions `beta + delta, beta * delta`. In fact, they can be just computed using
`find_identifiable_functions` function as we will explain in the next tutorial. Stay tuned!

[^1]: > A. Sedoglavic, [*A probabilistic algorithm to test local algebraic observability in polynomial time*](https://doi.org/10.1006/jsco.2002.0532), Journal of Symbolic Computation, 2002.
[^2]: > A. Ovchinnikov, A. Pillay, G. Pogudin, T. Scanlon, [*Multi-experiment Parameter Identifiability of ODEs and Model Theory*](https://doi.org/10.1137/21M1389845), SIAM Journal on Applied Algebra and Geometry, 2022.
[^3]: > D. Gonze, P. Ruoff, [*The Goodwin Oscillator and its Legacy*](https://doi.org/10.1007/s10441-020-09379-8), Acta Biotheoretica, 2020.
[^4]: > R. Dong, C. Goodbrake, H. Harrington, G. Pogudin, [*Differential elimination for dynamical models via projections with applications to structural identifiability*](https://doi.org/10.1137/22M1469067), SIAM Journal on Applied Algebra and Geometry, 2023.
