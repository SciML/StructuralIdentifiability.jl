## New release of SciML/StructuralIdentifiability.jl

On behalf of the package developers, I am glad to announce the 0.5.0 release of StructuralIdentifiability.jl, a toolbox for assessing the identifiability of ODEs.

### What is identifiability, and why do we care?

Identifiability is a property of ODEs and their parameters.

Inattention to identifiability can cause ambiguities in numerical software.

If a parameter is not identifiable, then it cannot be estimated from data, no matter how accurate the data is.
Thus, identifiability is a crucial property ... bla-bla-bla.

New release features:
- Finding all combinations of parameters that are identifiable  
- Automatic reparametrization of the model to make it identifiable

### Show me some cool code!

Let us walk through an example that involves a small real-world model.
In 2014, one possible reparametrization of this model was studied manually by Meshkat, Kuo, and DiStefano III in [^2]. 
Nowadays, in StructuralIdentifiability.jl, we can follow in their footsteps autonomosuly with a few lines of code.

In [^1], a study was conducted to measure the turnover of cholesterol in humans. A small amount of radioactive cholesterol was injected into the blood stream of a patient, and then the concentration of cholesterol in the blood stream and in the liver was measured over time.
Based on the data, the authors proposed a model, which, in a somewhat simplified form, is given by a system of ODEs:

$$
\left\{\begin{align*}
x_1'(t) &= k_{11} x_1 + k_{12} x_2 + u(t),\\
x_2'(t) &= k_{21} x_1 + k_{22} x_2,\\
y(t) &= x_1 / V
\end{align*}
\right.
$$

The physical interpretation of variables is as follows:
- The states $x_1(t)$ and $x_2(t)$ are the concentrations of cholesterol in the plasma and liver over time, respectively.
- The input $u(t)$ is the injection of cholesterol
- The output $y(t)$ is the quantity measured in the experiment, that is, the concentration of cholesterol in the plasma $x_1$ divided by the volume of plasma $V$.
- Some parameters $k_{11}, k_{12}, k_{21}, k_{22}$

To create this model in StructuralIdentifiability.jl, we simply write:

```julia
using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = k11 * x1 + k12 * x2 + u(t),
    x2'(t) = k21 * x1 + k22 * x2,
    y(t) = x1 / V
)
```

To assess identifiability of model variables, we can do

```julia
assess_identifiability(ode)
# returns
Dict{Any, Symbol} with 7 entries:
  k22 => :globally
  x2  => :nonidentifiable
  k21 => :nonidentifiable
  x1  => :globally
  k12 => :nonidentifiable
  k11 => :globally
  V   => :globally
```

Parameters ($k_{21}$ and $k_{12}$) are not identifiable, and the state $x_2$ is not identifiable.
Are there some combinations of parameters that are identifiable

```julia
find_identifiable_functions(ode, with_states = true)
# returns
6-element Vector{AbstractAlgebra.Generic.Frac{Nemo.fmpq_mpoly}}:
 x1
 V
 k11
 k22
 k12*x2
 k12*k21
```

So, although both $k_{21}$ and $k_{12}$ are not identifiable, their product $k_{12} k_{21}$ is identifiable.

```julia
new_ode, new_vars, _ = reparametrize_global(ode)

new_ode
# prints
X1'(t) = X1(t)*a2 + X2(t) + u1(t)
X2'(t) = X1(t)*a4 + X2(t)*a3
y1(t) = X1(t)//a1

new_vars
# prints
Dict{Nemo.fmpq_mpoly, AbstractAlgebra.Generic.Frac{Nemo.fmpq_mpoly}} with 8 entries:
  y1 => y
  u1 => u
  X1 => x1
  a4 => k12*k21
  a3 => k22
  X2 => k12*x2
  a1 => V
  a2 => k11
```

Important features of this reparametrization are:
- It is exact. No approximation errors were introduced.
- It makes everything identifiable
- It reduced the dimension: from 4 parameters to 3 parameters

<!--
In 2014, a similar reparametrization was obtained by hand by Meshkat, Kuo, and DiStefano III in [^2]. 
In this release of StructuralIdentifiability.jl, we can find a reparametrization of this model (and many other models) fully autonomosuly with a few lines of code.
-->

### Other stuff

Read more in our new preprint: [link](тык)

benchmarks[]: [link](тык)
documentation: [link](тык)

[^1]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC322312/
[^2]: https://pubmed.ncbi.nlm.nih.gov/25350289/
