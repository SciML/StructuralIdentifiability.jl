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
  title   = {Differential elimination for dynamical models via projections with applications to structural identifiability},
  journal = {preprint},
  url     = {https://arxiv.org/abs/2111.00991},
  year    = {2021}
}
```
## Feature Summary
`StructuralIdentifiability.jl` can assess local and global identifiability of ODE models. In addition to these straightforward identifiability queries on individual parameters, the package is able to distinguish between single- and multi-experiment identifiability.
## Feature List
* Local identifiability checks
* Global identifiability checks
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
Pkg.status(;mode = PKGMODE_MANIFEST) # hide
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
version = TOML.parse(read("../../Project.toml",String))["version"]
name = TOML.parse(read("../../Project.toml",String))["name"]
link = "https://github.com/SciML/"*name*".jl/tree/gh-pages/v"*version*"/assets/Manifest.toml"
```
```@raw html
">manifest</a> file and the
<a href="
```
```@eval
using TOML
version = TOML.parse(read("../../Project.toml",String))["version"]
name = TOML.parse(read("../../Project.toml",String))["name"]
link = "https://github.com/SciML/"*name*".jl/tree/gh-pages/v"*version*"/assets/Project.toml"
```
```@raw html
">project</a> file.
```