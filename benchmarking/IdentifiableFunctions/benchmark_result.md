## Benchmark results

2023-09-12T19:39:13.151

- Benchmarked function: `find_identifiable_functions`
- Workers: 8
- Timeout: 600 s

**All timings in seconds.**

|Model|with_states_VanDerHoevenLecerf / Runtime|with_states_VanDerHoevenLecerf / # Points, degree|with_states_VanDerHoevenLecerf / # Points, interpolation|with_states_CuytLee / Runtime|with_states_CuytLee / # Points, degree|with_states_CuytLee / # Points, interpolation|
|-----|---|---|---|---|---|---|
|MAPK model (5 outputs)|30.88|18|12|30.89|18|12|
|MAPK model (6 outputs)|5.53|18|12|5.52|18|12|
|Modified LV for testing|0.07|20|32|0.09|20|32|

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
