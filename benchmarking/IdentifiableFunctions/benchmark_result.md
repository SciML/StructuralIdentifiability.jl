## Benchmark results

2023-09-10T18:53:06.027

- Benchmarked function: `find_identifiable_functions`
- Workers: 8
- Timeout: 600 s

**All timings in seconds.**

|Model|(:gb,)_with_states / Runtime|(:gb,)_with_states / Polynomial?|(:hybrid, 1)_with_states / Runtime|(:hybrid, 1)_with_states / Polynomial?|(:hybrid, 2)_with_states / Runtime|(:hybrid, 2)_with_states / Polynomial?|(:hybrid, 3)_with_states / Runtime|(:hybrid, 3)_with_states / Polynomial?|
|-----|---|---|---|---|---|---|---|---|
|CGV1990|2.09|no|6.10|no| - | - |10.97|no|
|Goodwin oscillator|0.40|no|1.74|no|2.82|no|4.24|no|
|HIV2_io|1.95|no|14.76|no| - | - |35.25|no|
|KD1999|0.40|no|2.36|no|3.56|no|5.09|no|
|LLW1987_io|0.21|no|0.92|no|1.51|no|1.93|no|
|SEAIJRC Covid model|51.84|no| - | - | - | - | - | - |
|SEUIR|0.25|no|1.16|no| - | - |2.47|no|
|SLIQR|1.20|no|4.45|no| - | - |21.80|no|

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
* AbstractAlgebra : 0.27.10
* MacroTools : 0.5.11
* Nemo : 0.32.7
* SpecialFunctions : 2.3.1
