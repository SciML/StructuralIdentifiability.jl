# StructuralIdentifiability.jl

`StructuralIdentifiability.jl` is a comprehensive toolbox for assessing identifiability parameters.

This documentation contains information about the functionality of the package as well as examples of use cases.

## Installation

To install StructuralIdentifiability.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("StructuralIdentifiability")
```

## Citation
```latex
@article{structidjl,
  author  = {Dong, R. and Goodbrake, C. and Harrington, H. and Pogudin G.},
  title   = {Structural identifiability via input-output projections},
  journal = {Manuscript in preparation},
  year    = {2021}
}
```
## Feature Summary
`StructuralIdentifiability.jl` can assess local and global identifiability of ODE models. In addition to these straightforward identifiability queries on individual parameters, the package is able to distinguish between single- and multi-experiment identifiability.
## Feature List
* Local identifiability checks (single- and multi-experiment)
* Global identifiability checks (single- and multi-experiment)
* Assessment of identifiable functions of parameters
## Contributing

- Please refer to the
  [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
  for guidance on PRs, issues, and other matters relating to contributing to StructuralIdentifiability.
- There are a few community forums:
    - The #diffeq-bridged channel in the [Julia Slack](https://julialang.org/slack/)
    - [JuliaDiffEq](https://gitter.im/JuliaDiffEq/Lobby) on Gitter
    - On the Julia Discourse forums
    - See also [SciML Community page](https://sciml.ai/community/)
