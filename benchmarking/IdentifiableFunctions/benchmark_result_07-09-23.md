## Benchmark results

Timestamp: 2023-09-07T10:22:29.441

Timeout: 600 s

**Timings in seconds.**

|Model|Total time|Algebraic relations|Dim. before|Dim. after|
|-----|---|---|---|---|
|Akt pathway|92.44| no | 25 | 20 |
|Bilirubin2_io| - | - | - | - |
|Biohydrogenation_io|90.87| no | 10 | 9 |
|Bruno2016|90.93| no | 13 | 6 |
|CD8 T cell differentiation|90.68| no | 18 | 17 |
|CGV1990| - | - | - | - |
|Chemical reaction network|92.15| no | 12 | 12 |
|Crauste_SI|92.68| no | 18 | 17 |
|Fujita|99.71| no | 25 | 20 |
|Goodwin oscillator|94.11| no | 11 | 9 |
|HIV|96.67| no | 15 | 13 |
|HIV2_io| - | - | - | - |
|HighDimNonLin|119.25| no | 42 | 42 |
|JAK-STAT 1| - | - | - | - |
|KD1999|97.16| no | 19 | 14 |
|LLW1987_io|110.63| yes | 7 | 7 |
|LeukaemiaLeon2021| - | - | - | - |
|MAPK model (5 outputs bis)| - | - | - | - |
|MAPK model (5 outputs)|194.02| no | 34 | 34 |
|MAPK model (6 outputs)|135.00| no | 34 | 34 |
|Modified LV for testing|121.90| no | 6 | 5 |
|PK1| - | - | - | - |
|PK2|252.67| no | 11 | 11 |
|Pharm|259.99| no | 11 | 11 |
|Phosphorylation|125.16| no | 12 | 12 |
|Pivastatin|133.78| no | 11 | 10 |
|QWWC| - | - | - | - |
|QY| - | - | - | - |
|Ruminal lipolysis|90.83| no | 8 | 8 |
|SEAIJRC Covid model|192.52| no | 14 | 13 |
|SEIR 34|93.99| no | 12 | 10 |
|SEIR 36 ref|96.31| no | 21 | 20 |
|SEIR2T|90.70| no | 8 | 8 |
|SEIRT|95.07| no | 8 | 7 |
|SEIR_1_io|91.16| no | 9 | 8 |
|SEUIR|90.40| no | 10 | 7 |
|SIR 19|92.89| no | 11 | 10 |
|SIR 21|92.01| no | 11 | 10 |
|SIR 24|90.36| no | 9 | 6 |
|SIR 6|90.56| no | 7 | 5 |
|SIRS forced|101.59| no | 11 | 10 |
|SIWR original|103.76| no | 11 | 11 |
|SIWR with extra output|91.83| no | 11 | 11 |
|SLIQR| - | - | - | - |
|St| - | - | - | - |
|Transfection_4State|91.17| no | 9 | 8 |
|Treatment_io|92.82| no | 9 | 8 |
|TumorHu2019| - | - | - | - |
|TumorPillis2007| - | - | - | - |
|cLV1 (1o)| - | - | - | - |
|cLV1 (2o)|96.89| no | 18 | 17 |
|cholera|91.38| no | 11 | 11 |
|generalizedLoktaVolterra (1o)|93.30| no | 8 | 7 |
|generalizedLoktaVolterra (2o)|88.99| no | 8 | 8 |
|p53|91.15| no | 27 | 27 |

*Benchmarking environment:*

* Total RAM (GiB): 188
* Processor: Intel(R) Xeon(R) Gold 6130 CPU @ 2.10GHz
* Julia version: 1.9.2

Versions of the dependencies:

* Primes : 0.5.4
* BenchmarkTools : 1.3.2
* IterTools : 1.8.0
* PrecompileTools : 1.1.2
* Symbolics : 5.5.3
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
