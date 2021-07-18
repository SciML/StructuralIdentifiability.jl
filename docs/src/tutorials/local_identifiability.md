# Local Identifiability of Differential Models

In this tutorial, we will go over an example of solving a local identifiability problem for a simple system of ordinary differential equations.

We will introduce how to use the input parsing in `StructuralIdentifiability.jl` and the local identifiability assessment functionality.

## Input System

We will consider a simple two-species competition model

$x'_1 = k \,(1 - x_1 - x_2)\\ x'_2=r\,(1-x_1-x_2).$

To make it a proper input for the algorithm, we add an output function $y=x_1$ that equals to the population density of species 1 at any time $t$.

### Using the `@ODEmodel` macro

To parse the system of ordinary differential equations as above, we will use `@ODEmodel` macro. This is the easiest way to do so.

We have two state variables `x1, x2` (population densities), two parameters `k, r` (intrinsic growth rates), and one output function `y`. Note that there must be `(t)` to indicate time-dependent functions.

After using the macro, we use `assess_local_identifiability` function for that. This function accepts the ODE model, the probability of correctness, and the type of identifiability we would like to inquire about.

```@example
using StructuralIdentifiability

ode = @ODEmodel(
	x1'(t) = k * (1 - x1(t) - x2(t)),
	x2'(t) = r * (1 - x1(t) - x2(t)),
	y(t) = x1(t)
)

local_id = assess_local_identifiability(ode, 0.99)
```

The result shows that only the state variable's initial value $x'_1(0)$ is locally identifiable.

Let us now add another output function `y2(t)`:
```@example
using StructuralIdentifiability

ode = @ODEmodel(
	x1'(t) = k * (1 - x1(t) - x2(t)),
	x2'(t) = r * (1 - x1(t) - x2(t)),
	y1(t) = x1(t),
	y2(t) = x2(t) # new output function!
)

local_id = assess_local_identifiability(ode, 0.99) # this is a different result!
```

As you can see, for this new model with an additional output, all parameters are reported as locally identifiable with probability 0.99. 

## Note on Probability of Correctness
We set the probability of correctness $p$ to be `0.99`. Why would we ever want a lower value? Great question! The underlying algorithm relies on operations being modulo a large enough prime characteristic $\mathcal{P}\geq \kappa p$ where $\kappa$ is determined by the algorithm internally.

The algorithm's complexity is proportional to the size of operands (see proposition 3.1 in [[1]](#local)) and hence high probability of correctness may lead to higher size of coefficients during computation for some systems hence one may wish to lower $p$ to save on runtime (though in practice this is _very_ rare).

## References

<a id="local">[1]</a> A. Sedoglavic,
[*A probabilistic algorithm to test local algebraic observability in polynomial time*](https://doi.org/10.1006/jsco.2002.0532),
Journal of Symbolic Computation, 2002.

