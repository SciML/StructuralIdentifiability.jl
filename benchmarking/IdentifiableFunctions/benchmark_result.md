## Benchmark results

2023-09-11T16:22:23.532

- Benchmarked function: `find_identifiable_functions`
- Workers: 8
- Timeout: 600 s

**All timings in seconds.**

|Model|VanDerHoevenLecerf / Runtime|VanDerHoevenLecerf / # Points, degree|VanDerHoevenLecerf / # Points, interpolation|CuytLee / Runtime|CuytLee / # Points, degree|CuytLee / # Points, interpolation|
|-----|---|---|---|---|---|---|
|Bilirubin2_io|3.62|56|1656|4.47|56|1656|
|Biohydrogenation_io|0.38|72|624|0.50|72|624|
|CD8 T cell differentiation|0.30|20|20|0.37|20|20|
|Crauste_SI|0.36|20|20|0.25|20|20|
|HIV2_io|1.01|50|592|1.22|50|592|
|QY|30.66|156|7672|23.47|156|3836|
|SLIQR|1.39|76|1072|1.10|76|536|
|St|27.77|56|1472|28.98|56|1472|
|Treatment_io|0.48|72|264|0.55|72|264|

*Benchmarking environment:*

* Total RAM (GiB): 15
* Processor: Intel Xeon Processor (Icelake)
* Julia version: 1.9.3

Versions of the dependencies:

* Primes : 0.5.4
* BenchmarkTools : 1.3.2
* IterTools : 1.8.0
* PrecompileTools : 1.2.0
* Symbolics : 5.5.3
* Combinatorics : 1.0.2
* SymbolicUtils : 1.3.0
* DataStructures : 0.18.15
* Groebner : 0.4.3
* ParamPunPam : 0.0.3
* ModelingToolkit : 8.68.0
* AbstractAlgebra : 0.31.1
* MacroTools : 0.5.11
* Nemo : 0.35.3
* SpecialFunctions : 2.3.1
