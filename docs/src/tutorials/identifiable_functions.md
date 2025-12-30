# Globally Identifiable Functions

In addition to assessing identifiability of given functions of parameters or states, `StructuralIdentifiability.jl`
provides the function `find_identifiable_functions` which finds a set of identifiable functions such that any other
identifiable function can be expressed using them.
This allows to find out what actually is identifiable and what non-identifiability in the model at hand looks like.

For example, consider the following model[^1].

```@example funcs
using StructuralIdentifiability # hide
LLW1987 = @ODEmodel(
    x1'(t) = -p1 * x1(t) + p2 * u(t),
    x2'(t) = -p3 * x2(t) + p4 * u(t),
    x3'(t) = -(p1 + p3) * x3(t) + (p4 * x1(t) + p2 * x2(t)) * u(t),
    y1(t) = x3(t)
)
```

Several decades ago, this model was introduced to demonstrate the presence of nontrivial **un**identifiability in nonlinear systems of ODEs.
Nowadays, we can automatically find the identifiable combinations of parameters:

```@example funcs
using StructuralIdentifiability # hide
find_identifiable_functions(LLW1987)
```

From these expressions, we see that the values of `p1` and `p3` are not identifiable but an unordered pair
of numbers `{p1, p3}` is uniquely determined since `p1 + p3` and `p1 * p3` are known.
Furthermore, we see that, for fixed input and output, `p2` and `p4` can take infinitely many values but
knowing one of them, we would also be able to determine the other.

Moreover, we can find generators of all identifiable functions in parameters and states:

```@example funcs
find_identifiable_functions(LLW1987, with_states = true)
```

By default, `find_identifiable_functions` tries to simplify the output functions as much as possible, and it has `simplify` keyword responsible for
the degree of simplification. The default value is `:standard` but one could use `:strong` to try to simplify further
(at the expense of heavier computation) or use `:weak` to simplify less (but compute faster).

As `assess_identifiability` and `assess_local_identifiability`, `find_identifiable_functions` accepts an optional parameter `loglevel` (default: `Logging.Info`)
to adjust the verbosity of logging.

Finally, as for `assess_identifiability`, an experimental feature allows to provide an additional keyword argument `known_ic` to indicate functions of states and parameters for which the
initial conditions are assumed to be known (while the initial conditions of the system are still assumed to be generic).
In this case, the function will find identifiable functions of parameters and initial conditions rather than of parameters and states.
Let us add an assumption that the initial condition `x1(0)` is known:

```@example funcs
find_identifiable_functions(LLW1987, known_ic = [x1])
```

We see that `x2(0)` becomes an identifiable function as well, which is natural since `x1(t) * x2(t)` was an identifiable function before.

[^1]: > Y. Lecourtier, F. Lamnabhi-Lagarrigue, and E. Walter, [*A method to prove that nonlinear models can be unidentifiable*](https://doi.org/10.1109/CDC.1987.272467), Proceedings of the 26th Conference on Decision and Control, December 1987, 2144-2145;
