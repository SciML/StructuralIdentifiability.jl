## Benchmark results

2023-08-07T00:33:52.021


|Model|(:gb,)|(:gb,)_with_states|(:normalforms, 3)|(:normalforms, 3)_with_states|(:hybrid,)_with_states|(:hybrid,)|
|-----|---|---|---|---|---|---|
|Akt pathway| - | - | - | - | - | - |
|Bilirubin2_io|1.83|6.04|3.79|64.95| - |5.66|
|Biohydrogenation_io|13.78|5.35|15.60|7.08|4.81|11.90|
|CD8 T cell differentiation|5.54|11.90|7.73|25.40|21.48|6.14|
|Chemical reaction network|13.92|6.22|10.07|6.58|7.04|12.43|
|Fujita|15.97|14.52|18.52|57.84| - |13.37|
|Goodwin oscillator|1.69|2.05|2.16|2.57|3.32|2.97|
|HIV|6.92|15.75|6.35|15.61|15.45|5.48|
|HIV2_io|5.71|18.93|5.27|17.70|14.92|6.42|
|JAK-STAT 1|27.75|27.65|37.05| - | - |28.44|
|LLW1987_io|9.05|2.42|11.00|3.03|2.57|8.98|
|MAPK model (5 outputs bis)| - | - | - | - | - | - |
|MAPK model (5 outputs)|58.77| - |84.11| - | - |52.11|
|MAPK model (6 outputs)|19.46|24.55|46.46| - | - |39.66|
|Modified LV for testing|11.52|1.93|13.23|2.05|1.86|10.18|
|PK1|6.18|14.88|6.22|15.00|13.01|7.16|
|PK2|54.04| - |63.98| - | - |59.97|
|Pharm|61.19| - |45.15| - | - |45.61|
|QWWC|313.22| - |310.76| - | - |314.66|
|QY|9.58| - |9.16| - | - |8.37|
|SEAIJRC Covid model|66.22| - |80.82| - | - |52.03|
|SEIR 34|15.09| - |14.44| - | - |12.92|
|SEIR 36 ref|6.33| - |4.89| - | - |8.66|
|SEIR_1_io|9.00|2.94|8.46|4.45|3.65|11.66|
|SIR 19|14.33| - |10.38| - | - |13.66|
|SIR 21|17.40| - |16.29| - | - |17.65|
|SIR 24|8.70|2.74|7.52|3.73|3.96|8.20|
|SIR 6|13.13| - |12.33| - | - |12.10|
|SIRS forced|18.85|6.98|16.11|7.27|8.36|15.18|
|SIWR original|12.70|16.07|14.90|12.81|17.34|11.64|
|SIWR with extra output|5.58|7.05|4.75|5.99|7.03|4.42|
|SLIQR|9.00|8.89|12.51|12.48|12.73|12.13|
|St|34.64| - |34.52| - | - | - |
|Treatment_io|10.67| - |14.56| - | - |9.78|
|TumorHu2019| - | - | - | - | - | - |
|TumorPillis2007| - | - | - | - | - | - |

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
