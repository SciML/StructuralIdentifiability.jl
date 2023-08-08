## Benchmark results

2023-08-08T17:53:05.631

Timeout: 6000 s
|Model|(:gb,)|(:gb,)_with_states|(:normalforms, 3)|(:normalforms, 3)_with_states|(:hybrid,)|(:hybrid,)_with_states|
|-----|---|---|---|---|---|---|
|Akt pathway|11.31|11.89|13.42|58.36|14.54|64.56|
|Bilirubin2_io|1.53|4.49|2.17|66.33|6.58|140.03|
|Biohydrogenation_io|3.84|3.51|3.86|4.16|4.54|4.85|
|CD8 T cell differentiation|3.79|10.23|5.26|16.65|5.50|18.95|
|Chemical reaction network|3.71|3.92|3.64|5.36|4.36|6.69|
|Fujita|11.19|11.77|13.10|55.50|14.11|61.30|
|Goodwin oscillator|1.05|1.21|1.24|1.81|1.65|3.07|
|HIV|4.07|10.32|4.18|13.16|5.04|14.43|
|HIV2_io|3.98|12.59|4.36|14.70|6.01|22.18|
|JAK-STAT 1|24.21|27.25|28.86|294.78|30.32|299.23|
|LLW1987_io|0.95|1.52|1.14|1.75|1.45|2.55|
|MAPK model (5 outputs bis)| - | - | - | - | - | - |
|MAPK model (5 outputs)|65.54| - |76.97| - |71.03| - |
|MAPK model (6 outputs)|16.29|20.15|28.44|199.34|33.35|241.64|
|Modified LV for testing|0.88|1.32|1.06|1.61|1.31|2.25|
|PK1|3.73|10.43|4.05|14.89|8.14|21.62|
|PK2|72.49| - |78.18| - |102.62| - |
|Pharm|122.89| - |120.41| - |102.65| - |
|QWWC|764.57|1323.03|764.56|1268.76|763.13|902.24|
|QY|12.19| - |16.88| - |18.82| - |
|SEAIJRC Covid model|80.31|141.94|117.69|118.38|94.20|115.98|
|SEIR 34|5.36|32.87|20.60|60.35|32.12|63.33|
|SEIR 36 ref|4.21|13.04|5.57|51.85|5.60|118.83|
|SEIR_1_io|1.14|1.95|1.30|2.91|1.96|4.17|
|SIR 19|4.09|4.35|4.65|5.99|5.09|7.08|
|SIR 21|3.97|4.25|4.83|6.13|5.04|6.46|
|SIR 24|1.69|1.95|1.87|2.56|2.50|3.28|
|SIR 6|3.81|4.16|3.83|4.58|4.55|5.08|
|SIRS forced|6.98|6.75|7.43|10.11|8.68|9.38|
|SIWR original|13.77|12.56|11.78|13.51|12.28|14.76|
|SIWR with extra output|3.70|4.27|4.54|5.86|4.68|6.81|
|SLIQR|1.77|5.11|2.90|12.68|3.17|16.30|
|St|32.39|326.97|33.89|332.53| - | - |
|Treatment_io|3.26|3.77|3.50|4.98|4.50|6.08|
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
