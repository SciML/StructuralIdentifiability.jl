## Benchmark results

2023-08-08T15:08:30.629

Timeout: 6000 s
|Model|(:gb,)|(:gb,)_with_states|(:normalforms, 3)|(:normalforms, 3)_with_states|(:hybrid,)|(:hybrid,)_with_states|
|-----|---|---|---|---|---|---|
|Akt pathway|11.12|11.64|13.04|58.90|14.51|66.82|
|Bilirubin2_io|1.49|4.11|2.13|70.19|7.34| - |
|Biohydrogenation_io|3.82|3.36|3.89|4.17|4.78|5.15|
|CD8 T cell differentiation|3.73|10.37|4.88|17.74|5.63|19.56|
|Chemical reaction network|3.71|3.94|3.83|5.13|4.17|6.22|
|Fujita|10.58|11.69|13.15|55.88|13.64|62.97|
|Goodwin oscillator|1.09|1.26|1.29|1.92|1.84|2.94|
|HIV|4.09|10.61|4.42|12.74|4.87|14.62|
|HIV2_io|3.87|12.90|4.41|14.79|6.34|19.20|
|JAK-STAT 1|25.04|25.55|32.58|297.07|34.08|299.32|
|LLW1987_io|0.91|1.50|1.04|1.81|1.37|2.98|
|MAPK model (5 outputs bis)| - | - | - | - | - | - |
|MAPK model (5 outputs)|74.36| - |84.71| - |81.73| - |
|MAPK model (6 outputs)|24.06|35.51|43.11|183.27|37.57|178.21|
|Modified LV for testing|0.92|1.37|1.03|1.62|1.37|4.40|
|PK1|4.06|15.77|8.56|33.21|11.23|40.74|
|PK2|134.20| - |132.29| - |132.27| - |
|Pharm|124.90| - |121.70| - |115.34| - |
|QWWC|790.42|1361.00|824.03|1340.29|809.50|894.43|
|QY|26.65| - |31.69| - |1697.95| - |
|SEAIJRC Covid model|175.64| - |152.87| - |173.16| - |
|SEIR 34|4.97| - |6.63| - |17.38| - |
|SEIR 36 ref|4.72| - |7.61| - |10.42| - |
|SEIR_1_io|1.08|1.83|1.46|3.44|2.50|4.85|
|SIR 19|3.69| - |4.51| - |5.55| - |
|SIR 21|3.90| - |4.26| - |4.99| - |
|SIR 24|1.56|1.72|1.81|2.69|2.19|2.96|
|SIR 6|3.40| - |4.22| - |4.25| - |
|SIRS forced|6.03|6.83|7.05|8.90|7.39|9.45|
|SIWR original|11.18|12.74|11.33|13.38|12.53|14.13|
|SIWR with extra output|3.80|4.41|4.45|5.97|4.73|6.50|
|SLIQR|1.78|5.43|2.49|12.64|3.18|15.25|
|St|33.11|334.99|34.35|332.95| - | - |
|Treatment_io|3.26| - |3.50| - |4.51| - |
|TumorHu2019| - | - | - | - | - | - |
|TumorPillis2007| - | - | - | - | - | - |

*Benchmarking environment:*

* Total RAM (GiB): 188
* Processor: Intel(R) Xeon(R) Gold 6130 CPU @ 2.10GHz
* Julia version: 1.9.2

Versions of the dependencies:

* Primes : 0.5.4
* BenchmarkTools : 1.3.2
* IterTools : 1.8.0
* PrecompileTools : 1.1.2
* Symbolics : 5.5.1
* Combinatorics : 1.0.2
* SymbolicUtils : 1.2.0
* DataStructures : 0.18.15
* Groebner : 0.4.3
* ParamPunPam : 0.0.3
* ModelingToolkit : 8.64.0
* AbstractAlgebra : 0.27.10
* MacroTools : 0.5.10
* Nemo : 0.32.7
* SpecialFunctions : 2.3.0
