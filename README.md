# StructuralIdentifiability.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](http://si.sciml.ai/stable/)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/dev/modules/StructuralIdentifiability/)

[![codecov](https://codecov.io/gh/SciML/StructuralIdentifiability.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/StructuralIdentifiability.jl)
[![Build Status](https://github.com/SciML/StructuralIdentifiability.jl/workflows/CI/badge.svg)](https://github.com/SciML/StructuralIdentifiability.jl/actions?query=workflow%3ACI)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

## About

`StructuralIdentifiability.jl` is a Julia package for assessing structural parameter identifiability of parametric ODE models, both local and global.
For an introduction to structural identifiability, we refer to [[2]](#review).

## How to install

The package can be installed from this repository by
```julia
> using Pkg
> Pkg.add("StructuralIdentifiability")
```

## How to use

The package can be loaded by `using StructuralIdentifiability`.

### Defining a system

A parametric ODE system in the state-space from can be defined by the `@ODEmodel` macro:
```julia
ode = @ODEmodel(
    x1'(t) = -(a01 + a21) * x1(t) + a12 * x2(t) + u(t),
    x2'(t) = a21 * x1(t) - a12 * x2(t) - x3(t) / b,
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

For example, for the `ode` defined above, it will be
```julia
Dict{Nemo.fmpq_mpoly, Symbol} with 4 entries:
  a12 => :locally
  a21 => :globally
  a01 => :locally
  b   => :nonidentifiable
```

If one is interested in the identifiability of particular functions of the parameter, one can pass a list of them as a second argument:
```julia
assess_identifiability(ode, [a01 + a12, a01 * a12])
```
This will return a list of the results (i.e., `:globally`, `:locally`, or `:nonidentifiable`). In this example:
```julia
2-element Vector{Symbol}:
 :globally
 :globally
```

### Assessing local identifiability

Local identifiability can be assessed efficiently even for the models for which global identifiability analysis is out of reach. Moreover, the package can also assess local observability of the state variables. This can be done using the `assess_local_identifiability` function, for example:
```julia
assess_local_identifiability(ode)
```
The returned value is a dictionary from parameters and state variables to `1` (is locally identifiable/observable) and `0` (not identifiable/observable) values. In our example:
```julia
Dict{Nemo.fmpq_mpoly, Bool} with 7 entries:
  a12 => 1
  a21 => 1
  x3  => 0
  a01 => 1
  x2  => 1
  x1  => 1
  b   => 0
```

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

## Contacts

Maintained by Gleb Pogudin (<gleb.pogudin@polytechnique.edu>)

## References

<a id="global">[1]</a> 
Ruiwen Dong, Christian Goodbrake, Heather Harrington, and Gleb Pogudin,
[*Differential elimination for dynamical models via projections with applications to structural identifiability*](https://arxiv.org/abs/2111.00991),
preprint, 2021.

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

## Other identifiability software

The following software can be used to assess both local and global identifiability

* [SIAN](https://github.com/pogudingleb/SIAN) is written in Maple, there is a [Julia version](https://github.com/alexeyovchinnikov/SIAN-Julia). There is also a [webapp](https://maple.cloud/app/6509768948056064) with extended functionality.
* [DAISY](https://daisy.dei.unipd.it/) a package for the Reduce computer algebra system
* [COMBOS](http://biocyb1.cs.ucla.edu/combos/), a webapp.

Some benchmarking results for them can be found in [this paper](https://doi.org/10.1093/bioinformatics/bty1069).

The following software can be used to assess local identifiability

* [STRIKE-GOLDD](https://sites.google.com/site/strikegolddtoolbox/) in Matlab
* [ObservabilityTest](https://github.com/sedoglavic/ObservabilityTest/) in Maple
* [IdentifiabilityAnalysis](http://www.fcc.chalmers.se/software/other-software/identifiabilityanalysis/) in Mathematica

If your software is not listed here, sorry, we either forgot or did not know about it, feel free to contact Gleb Pogudin.
