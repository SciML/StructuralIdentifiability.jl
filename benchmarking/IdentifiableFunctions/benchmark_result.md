## Benchmark results

2023-09-15T02:10:37.348

- Benchmarked function: `find_identifiable_functions`
- Workers: 7
- Timeout: 600 s

**All timings in seconds.**

|Model|(:normalforms, 2)_with_states / Runtime|(:normalforms, 2)_with_states / beautifulization|
|-----|---|---|
|Akt pathway|2.44|0.28|
|Bilirubin2_io|4.08|0.76|
|Biohydrogenation_io|0.44|0.11|
|Bruno2016|0.19|0.07|
|CD8 T cell differentiation|0.72|0.27|
|CGV1990|2.25|0.31|
|Chemical reaction network|0.36|0.13|
|Crauste_SI|0.47|0.18|
|Fujita|2.50|0.30|
|Goodwin oscillator|0.29|0.09|
|HIV|0.78|0.19|
|HIV2_io|2.00|0.27|
|HighDimNonLin|20.55|0.68|
|KD1999|0.97|0.14|
|LLW1987_io|0.25|0.08|
|MAPK model (5 outputs)|48.45|0.37|
|MAPK model (6 outputs)|7.99|0.42|
|Modified LV for testing|0.07|0.01|
|PK1|0.74|0.18|
|PK2|130.96|0.07|
|Pharm|152.26|0.09|
|Phosphorylation|0.45|0.11|
|Pivastatin|8.44|0.04|
|QY|49.21|0.23|
|Ruminal lipolysis|0.21|0.07|
|SEAIJRC Covid model|85.13|0.11|
|SEIR 34|0.39|0.09|
|SEIR 36 ref|0.78|0.14|
|SEIR2T|0.23|0.09|
|SEIRT|0.36|0.08|
|SEIR_1_io|0.56|0.10|
|SEUIR| - | - |
|SIR 19| - | - |
|SIR 21|0.30|0.19|
|SIR 24|0.36|0.07|
|SIR 6|0.20|0.08|
|SIRS forced|10.62|0.40|
|SIWR original|13.72|0.07|
|SIWR with extra output|0.67|0.10|
|SLIQR|1.42|0.15|
|St|37.54|1.38|
|Transfection_4State|0.23|0.07|
|Treatment_io|0.44|0.14|
|cLV1 (2o)|0.80|0.20|
|cholera|0.66|0.09|
|generalizedLoktaVolterra (1o)|0.18|0.05|
|generalizedLoktaVolterra (2o)|0.09|0.07|
|p53|1.21|0.34|

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
