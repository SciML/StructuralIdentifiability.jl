## Benchmark results

2023-08-06T06:03:27.527


|Model|(:gb,)|(:normalforms, 3)|(:normalforms, 2)|(:hybrid,)|
|-----|---|---|---|---|
|Akt pathway|14.36|11.89|11.62|12.85|
|Bilirubin2_io|1.67|2.46|2.57|13.02|
|Biohydrogenation_io|10.91|8.76|10.36|9.76|
|CD8 T cell differentiation|4.36|7.23|4.42|5.61|
|Chemical reaction network|11.34|8.92|10.98|9.48|
|Fujita|10.77|13.85|9.25|12.79|
|Goodwin oscillator|1.20|1.50|1.45|1.62|
|HIV|4.53|5.90|6.19|3.94|
|HIV2_io|4.29|4.73|4.53|7.41|
|JAK-STAT 1|23.80|27.03|20.56|24.45|
|LLW1987_io|8.66|8.41|8.70|6.24|
|MAPK model (5 outputs bis)| - | - | - | - |
|MAPK model (5 outputs)|43.64|49.23|44.36|46.11|
|MAPK model (6 outputs)|17.21|18.10|17.44|19.49|
|Modified LV for testing|7.95|8.22|8.64|8.30|
|PK1|4.19|3.55|4.29|4.86|
|PK2|49.98|46.94|51.42|45.81|
|Pharm|49.68|47.79|47.40|49.05|
|QWWC|313.22|310.76|312.59|314.66|
|QY|9.75| - | - | - |
|SEAIJRC Covid model|52.32|52.25|52.35|49.02|
|SEIR 34|10.79|8.70|10.78|10.59|
|SEIR 36 ref|4.44|5.09|4.78|4.45|
|SEIR_1_io|7.01|8.87|10.91|8.21|
|SIR 19|12.83|9.57|10.53|8.80|
|SIR 21|10.46|9.23|11.07|9.70|
|SIR 24|7.66|8.95|6.71|7.09|
|SIR 6|11.49|10.32|11.09|9.13|
|SIRS forced|11.45|13.40|14.64|12.13|
|SIWR original|11.30|11.20|8.63|10.19|
|SIWR with extra output|4.45|4.34|3.65|4.60|
|SLIQR|9.81|9.29|9.90|9.48|
|St|29.29|27.43|23.87| - |
|Treatment_io|10.72|10.97|11.16|10.64|
|TumorHu2019| - | - | - | - |
|TumorPillis2007| - | - | - | - |

*Benchmarking environment:*

* Total RAM (GiB): 2015
* Processor: AMD EPYC 7702 64-Core Processor                
* Julia version: 1.9.1

Versions of the dependencies:

* Primes : 0.5.4
* BenchmarkTools : 1.3.2
* IterTools : 1.8.0
* PrecompileTools : 1.1.2
* Symbolics : 5.5.0
* Combinatorics : 1.0.2
* SymbolicUtils : 1.1.0
* DataStructures : 0.18.14
* Groebner : 0.4.3
* ParamPunPam : 0.0.3
* ModelingToolkit : 8.63.0
* AbstractAlgebra : 0.27.10
* MacroTools : 0.5.10
* Nemo : 0.32.7
* SpecialFunctions : 2.3.0
