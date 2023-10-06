# StructuralIdentifiability.jl

`StructuralIdentifiability.jl` is a comprehensive toolbox for assessing identifiability parameters.

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
  title   = {Differential Elimination for Dynamical Models via Projections with Applications to Structural Identifiability},
  journal = {SIAM Journal on Applied Algebra and Geometry},
  url     = {https://doi.org/10.1137/22M1469067},
  year    = {2023}
  volume  = {7},
  number  = {1},
  pages   = {194-235}
}
```

## Feature Summary

`StructuralIdentifiability.jl` can assess local and global identifiability of ODE models. In addition to these straightforward identifiability queries on individual parameters, the package can distinguish between single- and multi-experiment identifiability.

## Feature List

  - Local identifiability checks
  - Global identifiability checks
  - Assessment of identifiable functions of parameters and states
  - Model reparametrization (experimental)

## Contributing

  - Please refer to the
    [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
    for guidance on PRs, issues, and other matters relating to contributing to StructuralIdentifiability.

  - There are a few community forums:
    
      + The #diffeq-bridged channel in the [Julia Slack](https://julialang.org/slack/)
      + [JuliaDiffEq](https://gitter.im/JuliaDiffEq/Lobby) on Gitter
      + On the Julia Discourse forums
      + See also [SciML Community page](https://sciml.ai/community/)

## Reproducibility

```@raw html
<details><summary>The documentation of this SciML package was built using these direct dependencies,</summary>
```

```@example
using Pkg # hide
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>and using this machine and Julia version.</summary>
```

```@example
using InteractiveUtils # hide
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```

```@example
using Pkg # hide
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```

```@raw html
You can also download the 
<a href="
```

```@eval
using TOML
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link =
    "https://github.com/SciML/" *
    name *
    ".jl/tree/gh-pages/v" *
    version *
    "/assets/Manifest.toml"
nothing
```

```@raw html
">manifest</a> file and the
<a href="
```

```@eval
using TOML
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link =
    "https://github.com/SciML/" *
    name *
    ".jl/tree/gh-pages/v" *
    version *
    "/assets/Project.toml"
nothing
```

```@raw html
">project</a> file.
```
