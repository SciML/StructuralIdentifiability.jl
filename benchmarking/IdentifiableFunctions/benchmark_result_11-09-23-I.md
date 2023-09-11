## Benchmark results

2023-09-11T11:34:49.665

- Benchmarked function: `find_identifiable_functions`
- Workers: 8
- Timeout: 300 s

**All timings in seconds.**

|Model|VanDerHoevenLecerf / Runtime|VanDerHoevenLecerf / # Points, degree|VanDerHoevenLecerf / # Points, interpolation|CuytLee / Runtime|CuytLee / # Points, degree|CuytLee / # Points, interpolation|
|-----|---|---|---|---|---|---|
|Akt pathway|1.21|134|136| - | - | - |
|Bilirubin2_io|4.11|136|1680|7.85|136|1680|
|Biohydrogenation_io|0.59|216|1176|1.27|216|1176|
|Bruno2016|0.17|36|28| - | - | - |
|CD8 T cell differentiation|0.27|36|20|0.34|36|20|
|CGV1990|3.35|264|168|3.17|264|168|
|Chemical reaction network|0.37|36|26|0.41|36|16|
|Crauste_SI|0.43|36|20|0.35|36|20|
|Fujita|2.01|134|136| - | - | - |
|Goodwin oscillator|0.38|52|48|0.39|52|48|
|HIV|0.54|102|40|0.52|102|40|
|HIV2_io|2.58|200|592|2.64|200|592|
|HighDimNonLin|21.45|36|22|20.69|36|22|
|JAK-STAT 1| - | - | - | - | - | - |
|KD1999|0.62|84|24|0.65|84|24|
|LLW1987_io|0.24|52|40| - | - | - |
|LeukaemiaLeon2021| - | - | - | - | - | - |
|MAPK model (5 outputs bis)| - | - | - | - | - | - |
|MAPK model (5 outputs)|61.68|84|16|63.70|84|16|
|MAPK model (6 outputs)|12.32|84|16| - | - | - |
|Modified LV for testing|0.06|36|32| - | - | - |
|PK1|0.44|68|96| - | - | - |
|PK2| - | - | - | - | - | - |
|Pharm| - | - | - | - | - | - |
|Phosphorylation|0.47|36|26|0.52|36|16|
|Pivastatin|8.58|52|22|8.62|52|22|
|QWWC| - | - | - | - | - | - |
|QY| - | - | - | - | - | - |
|Ruminal lipolysis|0.23|52|16|0.24|52|16|
|SEAIJRC Covid model| - | - | - | - | - | - |
|SEIR 34|0.41|136|124|0.56|136|124|
|SEIR 36 ref|1.12|52|20|1.30|52|22|
|SEIR2T|0.21|52|16|0.27|52|18|
|SEIRT|0.28|86|120|0.44|86|128|
|SEIR_1_io|0.50|102|168|0.47|102|192|
|SEUIR|0.60|68|36|0.47|68|46|
|SIR 19|0.30|36|18|0.27|36|18|
|SIR 21|0.26|36|18|0.27|36|18|
|SIR 24|0.52|84|48|0.37|84|52|
|SIR 6|0.19|86|36|0.17|86|38|
|SIRS forced|10.91|36|20|11.10|36|20|
|SIWR original|16.19|52|16|16.29|52|18|
|SIWR with extra output| - | - | - |0.79|52|18|
|SLIQR| - | - | - |1.39|170|624|
|St| - | - | - | - | - | - |
|Transfection_4State| - | - | - |0.18|36|18|
|Treatment_io| - | - | - |0.56|136|296|
|TumorHu2019| - | - | - | - | - | - |
|TumorPillis2007| - | - | - | - | - | - |
|cLV1 (1o)| - | - | - | - | - | - |
|cLV1 (2o)|0.51|52|20|0.33|52|22|
|cholera|0.74|52|16|0.73|52|18|
|generalizedLoktaVolterra (1o)|0.13|36|26|0.18|36|26|
|generalizedLoktaVolterra (2o)|0.15|36|12|0.14|36|12|
|p53|1.17|200|84|1.26|200|74|

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
