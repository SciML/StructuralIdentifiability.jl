## Benchmark results

2023-09-10T05:14:42.659

- Benchmarked function: `find_identifiable_functions`
- Workers: 4
- Timeout: 200 s

**All timings in seconds.**

|Model|with_states / Runtime|with_states / Polynomial?|(:normalforms, 3)_with_states / Runtime|(:normalforms, 3)_with_states / Polynomial?|(:hybrid, 1)_with_states / Runtime|(:hybrid, 1)_with_states / Polynomial?|
|-----|---|---|---|---|---|---|
|Akt pathway| - | - | - | - | - | - |
|Bilirubin2_io| - | - | - | - | - | - |
|Biohydrogenation_io| - | - | - | - | - | - |
|Bruno2016| - | - | - | - | - | - |
|CD8 T cell differentiation| - | - | - | - | - | - |
|CGV1990| - | - | - | - | - | - |
|Chemical reaction network| - | - | - | - | - | - |
|Crauste_SI| - | - | - | - | - | - |
|Fujita| - | - | - | - | - | - |
|Goodwin oscillator| - | - | - | - | - | - |
|HIV| - | - | - | - | - | - |
|HIV2_io| - | - | - | - | - | - |
|HighDimNonLin| - | - | - | - | - | - |
|JAK-STAT 1| - | - | - | - | - | - |
|KD1999| - | - | - | - | - | - |
|LLW1987_io| - | - | - | - | - | - |
|LeukaemiaLeon2021| - | - | - | - | - | - |
|MAPK model (5 outputs bis)| - | - | - | - | - | - |
|MAPK model (5 outputs)| - | - | - | - | - | - |
|MAPK model (6 outputs)| - | - | - | - | - | - |
|Modified LV for testing| - | - | - | - | - | - |
|NFkB| - | - | - | - | - | - |
|PK1| - | - | - | - | - | - |
|PK2| - | - | - | - | - | - |
|Pharm| - | - | - | - | - | - |
|Phosphorylation| - | - | - | - | - | - |
|Pivastatin| - | - | - | - | - | - |
|QWWC| - | - | - | - | - | - |
|QY| - | - | - | - | - | - |
|Ruminal lipolysis| - | - | - | - | - | - |
|SEAIJRC Covid model| - | - | - | - | - | - |
|SEIR 34| - | - | - | - | - | - |
|SEIR 36 ref| - | - | - | - | - | - |
|SEIR2T| - | - | - | - | - | - |
|SEIRT| - | - | - | - | - | - |
|SEIR_1_io| - | - | - | - | - | - |
|SEUIR| - | - | - | - | - | - |
|SIR 19| - | - | - | - | - | - |
|SIR 21| - | - | - | - | - | - |
|SIR 24| - | - | - | - | - | - |
|SIR 6| - | - | - | - | - | - |
|SIRS forced| - | - | - | - | - | - |
|SIWR original| - | - | - | - | - | - |
|SIWR with extra output| - | - | - | - | - | - |
|SLIQR| - | - | - | - | - | - |
|St| - | - | - | - | - | - |
|Transfection_4State| - | - | - | - | - | - |
|Treatment_io| - | - | - | - | - | - |
|TumorHu2019| - | - | - | - | - | - |
|TumorPillis2007| - | - | - | - | - | - |
|cLV1 (1o)| - | - | - | - | - | - |
|cLV1 (2o)| - | - | - | - | - | - |
|cholera| - | - | - | - | - | - |
|generalizedLoktaVolterra (1o)| - | - | - | - | - | - |
|generalizedLoktaVolterra (2o)| - | - | - | - | - | - |
|p53| - | - | - | - | - | - |

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
