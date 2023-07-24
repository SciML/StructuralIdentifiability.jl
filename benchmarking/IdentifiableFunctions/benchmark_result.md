## Benchmark results

2023-07-22T10:39:24.453


|Model|io|global id.|inclusion|inclusion Zp|ParamPunPam.jl|total|
|-----|---|---|---|---|---|---|
|Akt pathway|2.68|0.00|1.43|0.91|3.93|9.56|
|Bilirubin2_io|0.01|0.00|0.02|0.01|0.69|1.43|
|Biohydrogenation_io|2.82|0.00|1.38|0.92|3.38|9.08|
|CD8 T cell differentiation|2.56|0.00|0.00|0.00|0.68|3.80|
|Chemical reaction network|2.61|0.00|1.45|0.93|3.29|8.88|
|Fujita|2.57|0.00|1.37|0.88|3.77|9.18|
|Goodwin oscillator|0.01|0.00|0.03|0.00|0.30|0.95|
|HIV|2.79|0.00|0.01|0.04|0.79|4.30|
|HIV2_io|2.38|0.00|0.00|0.00|0.81|3.82|
|JAK-STAT 1|10.18|0.00|1.64|1.16|4.31|17.98|
|LLW1987_io|0.00|0.00|1.37|0.89|3.09|5.89|
|MAPK model (5 outputs bis)| - | - | - | - | - | - |
|MAPK model (5 outputs)|30.51|0.00|2.81|2.02|9.14|45.32|
|MAPK model (6 outputs)|5.41|0.00|1.72|1.18|4.46|13.45|
|Modified LV for testing|0.00|0.00|1.43|0.92|3.24|6.17|
|NFkB| - | - | - | - | - | - |
|PK1|2.40|0.00|0.00|0.00|0.72|3.72|
|PK2|14.13|0.00|18.66|4.55|10.60|50.83|
|Pharm|14.96|0.00|17.95|4.55|10.47|50.79|
|QWWC|156.45|0.00|119.86|26.51|39.48|363.14|
|QY|0.12|0.00|4.04|0.13|2.36|7.33|
|SEAIJRC Covid model|23.44|0.00|15.28|3.01|11.20|55.66|
|SEIR 34|2.47|0.00|1.38|0.89|3.19|8.54|
|SEIR 36 ref|2.67|0.00|0.01|0.01|0.74|4.04|
|SEIR_1_io|0.01|0.00|1.37|0.89|3.15|5.97|
|SIR 19|2.46|0.00|1.36|0.89|3.04|8.31|
|SIR 21|2.47|0.00|1.39|0.86|3.06|8.33|
|SIR 24|0.14|0.00|1.39|0.98|3.23|6.40|
|SIR 6|2.58|0.00|1.38|0.88|3.27|8.66|
|SIRS forced|3.87|0.00|1.82|1.22|3.89|11.68|
|SIWR original|2.56|0.00|3.23|0.77|1.88|9.40|
|SIWR with extra output|3.30|0.00|0.10|0.02|0.45|4.53|
|SLIQR|0.01|0.00|1.46|0.98|3.69|6.78|
|St|4.11|0.00|3.20|0.68|15.59|25.65|
|Treatment_io|2.49|0.00|1.45|0.94|3.32|8.80|
|TumorHu2019| - | - | - | - | - | - |
|TumorPillis2007| - | - | - | - | - | - |

*Benchmarking environment:*

* Total RAM (GiB): 2015
* Processor: AMD EPYC 7702 64-Core Processor                
* Julia version: 1.9.1

Versions of the dependencies:

* Primes : 0.5.3
* BenchmarkTools : 1.3.2
* IterTools : 1.8.0
* PrecompileTools : 1.1.2
* Symbolics : 5.5.0
* SymbolicUtils : 1.0.5
* DataStructures : 0.18.14
* Groebner : 0.3.6
* ParamPunPam : 0.0.1
* ModelingToolkit : 8.60.0
* AbstractAlgebra : 0.27.10
* MacroTools : 0.5.10
* Nemo : 0.32.7
* SpecialFunctions : 2.3.0
