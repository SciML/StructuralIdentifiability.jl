## Benchmark results

2023-09-10T01:08:41.878

Function: `find_identifiable_functions`

Workers: 4

Timeout: 100 s

**All timings in seconds.**

|Model|Runtime|
|-----|---|
|Akt pathway|77.06|
|Bilirubin2_io| - |
|Biohydrogenation_io| - |

*Benchmarking environment:*

* Total RAM (GiB): 7
* Processor: Intel(R) Core(TM) i5-8250U CPU @ 1.60GHz
* Julia version: 1.9.1

Versions of the dependencies:

* Primes : 0.5.3
* BenchmarkTools : 1.3.2
* IterTools : 1.8.0
* PrecompileTools : 1.1.2
* Symbolics : 5.5.1
* Combinatorics : 1.0.2
* SymbolicUtils : 1.1.0
* DataStructures : 0.18.14
* Groebner : 0.4.3
* ParamPunPam : 0.0.3
* ModelingToolkit : 8.62.0
* AbstractAlgebra : 0.27.10
* MacroTools : 0.5.10
* Nemo : 0.32.7
* SpecialFunctions : 2.3.0
