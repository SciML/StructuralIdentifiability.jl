## Benchmark results

2023-07-13T15:36:23.490


|Model|io|global id.|gen. ideal|inclusion|ParamPunPam.jl|total|
|-----|---|---|---|---|---|---|
|Akt pathway|3.26|4.66|0.03|3.84|4.65|16.66|
|Bilirubin2_io| - | - | - | - | - | - |
|Biohydrogenation_io|2.64|0.01|0.02|0.01|2.98|5.85|
|CD8 T cell differentiation|2.66|3.30|0.02|3.01|1.49|10.60|
|Chemical reaction network|2.66|0.02|0.02|0.05|2.72|5.93|
|Fujita|2.67|3.93|0.02|3.57|4.32|14.68|
|Goodwin oscillator|0.05|0.01|0.02|0.05|1.14|1.44|
|HIV|2.66|0.03|0.02|0.02|1.61|7.77|
|HIV2_io|2.47|0.01|0.02|0.01|1.58|7.61|
|LLW1987_io|0.00|2.90|0.02|0.00|2.93|6.02|
|MAPK model (5 outputs bis)| - | - | - | - | - | - |
|MAPK model (5 outputs)|32.16|18.35|0.18|156.47|225.27|433.32|
|MAPK model (6 outputs)|5.76|5.01|0.06|33.91|52.03|97.56|
|Modified LV for testing|0.00|3.14|0.02|0.00|3.32|6.69|
|NFkB| - | - | - | - | - | - |
|Pharm|14.76|16.96|2.05|499.75|119.09|655.68|
|QWWC| - | - | - | - | - | - |
|QY|0.12|0.44|0.02|5.39|5.64|12.09|
|SEAIJRC Covid model|23.28|13.25|0.99|34.85|33.59|108.83|
|SEIR_1_io|0.01|2.92|0.02|0.04|2.97|6.15|
|SIRS forced|4.05|0.61|0.07|1.47|4.24|10.69|
|SIWR original|2.64|3.08|0.28|62.59|271.53|341.05|
|SIWR with extra output|3.05|0.08|0.03|3.72|1.85|9.21|
|SLIQR|0.01|3.13|0.02|3.57|3.31|10.28|
|St| - | - | - | - | - | - |
|Treatment_io|2.46|2.86|0.02|0.05|3.01|8.54|

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
