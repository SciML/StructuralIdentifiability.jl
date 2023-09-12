## Benchmark results

2023-09-12T21:39:54.563

- Benchmarked function: `find_identifiable_functions`
- Workers: 8
- Timeout: 600 s

**All timings in seconds.**

|Model|with_states_VanDerHoevenLecerf / Runtime|with_states_VanDerHoevenLecerf / # Points, degree|with_states_VanDerHoevenLecerf / # Points, interpolation|with_states_CuytLee / Runtime|with_states_CuytLee / # Points, degree|with_states_CuytLee / # Points, interpolation|
|-----|---|---|---|---|---|---|
|Akt pathway|1.03|32|120|2.10|32|120|
|Bilirubin2_io|3.93|56|1656|8.12|56|1656|
|Biohydrogenation_io|0.64|72|624|0.97|72|624|
|Bruno2016|0.25|18|28|0.26|18|28|
|CD8 T cell differentiation|0.37|20|20|0.40|20|20|
|CGV1990|2.32|56|148|2.44|56|148|
|Chemical reaction network|0.38|18|26|0.37|18|26|
|Crauste_SI|0.38|20|20|0.40|20|20|
|Fujita|0.99|32|120|1.84|32|120|
|Goodwin oscillator|0.24|20|44|0.37|20|44|
|HIV|0.46|38|40|0.52|38|40|
|HIV2_io|1.97|50|592|2.22|50|592|
|HighDimNonLin|16.88|16|22|21.07|16|22|
|JAK-STAT 1| - | - | - | - | - | - |
|KD1999|0.55|20|22|0.65|20|22|
|LLW1987_io|0.20|20|40|0.19|20|40|
|LeukaemiaLeon2021| - | - | - | - | - | - |
|MAPK model (5 outputs bis)| - | - | - | - | - | - |
|MAPK model (5 outputs)|30.88|18|12|30.89|18|12|
|MAPK model (6 outputs)|5.53|18|12|5.52|18|12|
|Modified LV for testing|0.07|20|32|0.09|20|32|
|PK1|0.28|20|96|0.40|20|96|
|PK2|201.37|20|12|204.41|20|12|
|Pharm|188.66|20|12|201.56|20|12|
|Phosphorylation|0.44|18|26|0.48|18|26|
|Pivastatin|8.17|20|20|8.56|20|20|
|QWWC| - | - | - | - | - | - |
|QY|30.66|156|7672|23.47|156|3836|
|Ruminal lipolysis|0.13|18|12|0.20|18|12|
|SEAIJRC Covid model|51.42|20|44|87.08|20|22|
|SEIR 34|0.26|56|120|0.47|56|120|
|SEIR 36 ref|0.70|20|16|1.05|20|16|
|SEIR2T|0.14|18|12|0.21|18|12|
|SEIRT|0.29|32|120|0.30|32|120|
|SEIR_1_io|0.40|38|168|0.42|38|168|
|SEUIR|0.32|20|34|0.35|20|34|
|SIR 19|0.26|18|18|0.27|18|18|
|SIR 21|0.24|18|18|0.26|18|18|
|SIR 24|0.29|20|40|0.33|20|40|
|SIR 6|0.14|36|34|0.14|36|34|
|SIRS forced|10.88|20|20|10.59|20|20|
|SIWR original|14.05|18|12|14.82|18|12|
|SIWR with extra output|0.65|18|12|0.63|18|12|
|SLIQR|1.39|76|1072|1.10|76|536|
|St|27.77|56|1472|28.98|56|1472|
|Transfection_4State|0.12|18|18|0.20|18|18|
|Treatment_io|0.48|72|264|0.55|72|264|
|TumorHu2019| - | - | - | - | - | - |
|TumorPillis2007| - | - | - | - | - | - |
|cLV1 (1o)| - | - | - | - | - | - |
|cLV1 (2o)|0.35|18|20|0.46|18|20|
|cholera|0.64|18|12|0.62|18|12|
|generalizedLoktaVolterra (1o)|0.14|18|26|0.15|18|26|
|generalizedLoktaVolterra (2o)|0.16|16|12|0.16|16|12|
|p53|0.96|56|64|1.01|56|44|

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
