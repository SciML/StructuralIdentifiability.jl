# StructuralIdentifiability.jl

## About

`StructuralIdentifiability.jl` is a Julia package for assessing structural parameter identifiability of parametric ODE models, both local and global.
For an introduction to structural identifiability, we refer to [[2]](#review).

## How to install

The package can be installed from this repository by
```julia
> using Pkg
> Pkg.add("https://github.com/pogudingleb/StructuralIdentifiability.jl")
```

## How to use

The package can be loaded by `using StructuralIdentifiability`.

### Defining a system

A parametric ODE system in the state-space from can be defined by the `@ODEmodel` macro:
```julia
ode = @ODEmodel(
    x1'(t) = -(a01 + a21) * x1(t) + a12 * x2(t) + u(t),
    x2'(t) = a21 * x1(t) - a21 * x2(t) - x3(t) / b,
    x3'(t) = x3(t),
    y(t) = x2(t)
)
```
In this example:

* `x1(t), x2(t), x3(t)` are the **state variables**, they defined the state of the system and are assumed to be unknown;
* `u(t)` is the **input/control variable** which is assumed to be known and generic (exciting) enough;
* `y(t)` is the **output variable** which is assumed to be observed in the experiments and, thus, known;
* `a01, a21, a12, b` are unknown scalar **parameters**.

Note that there may be mulitple inputs and outputs.

### Assessing identifiability

The identifiability of the parameters in the model can be assessed by the `assess_identifiability` function as follows
```julia
assess_identifiability(ode)
```
The returned value is a dictionary from the parameter of the model to one of the symbols 

* `:globally` meaning that the parameter is globally identifiable
* `:locally` meaning that the parameter is locally but not globally identifiable
* `:nonidentifiable` meaning that the parameter is not identifiable even locally.

If one is interested in the identifiability of particular functions of the parameter, one can pass a list of them as a second argument:
```julia
assess_identifiability(ode, [a01 + a21, a01 * a21])
```
This will return a list of the results (i.e., `:globally`, `:locally`, or `:nonidentifiable`).

### Assessing local identifiability

Local identifiability can be assessed efficiently even for the models for which global identifiability analysis is out of reach. Moreover, the package can also assess local observability of the state variables. This can be done using the `assess_local_identifiability` function, for example:
```julia
assess_local_identifiability(ode)
```
The returned value is a dictionary from parameters and state variables to `true` (is locally identifiable) and `false` (not identifiable) values.
As for `assess_identifiability`, one can assess local identifiability of arbitrary rational functions in the parameters (and also states) by providing a list of such functions as the second argument.

**Remark** The algorithms we used are randomized, the default probability of the correctness of the result is 99%, one can change it by changing the value of a keyword argument `p` to any real number between 0 and 1, for example:
```julia
# pobability of correctness 99.9%
assess_identifiability(ode; p=0.999)
```

## Algorithms

The algorithm used for assessing global identifiability is described in [[1]](#global). 
Local identifiability is assessed using the algorithm by Sedoglavic [[4]](#local).
We also use some of the algorithms described in [[3]](#allident).

## References

<a id="global">[1]</a> 
Ruiwen Dong, Christian Goodbrake, Heather Harrington, and Gleb Pogudin,
*Structural identifiability via input-output projections*,
in preparation.

<a id="review">[2]</a> 
Hongyu Miao, Xiaohua Xia, Alan S. Perelson, and Hulin Wu,
[*On Identifiability of Nonlinear ODE Models and Applications in Viral Dynamics*](https://doi.org/10.1137/090757009),
SIAM Review, 2011.

<a id="allident">[3]</a> 
Alexey Ovchinnikov, Anand Pillay, Gleb Pogudin, and Thomas Scanlon,
[*Computing all identifiable functions for ODE models*](https://arxiv.org/abs/2004.07774),
preprint, 2020.

<a id="local">[4]</a> 
Alexandre Sedoglavic,
[*A probabilistic algorithm to test local algebraic observability in polynomial time*](https://doi.org/10.1006/jsco.2002.0532),
Journal of Symbolic Computation, 2002.

