## Simplification scores

*Smaller is better.*

2023-08-10T10:26:34.628

Timeout: 2000 s
|Model|(:gb,)_with_states|(:normalforms, 2)_with_states|(:normalforms, 3)_with_states|(:gbfan, 10)_with_states|(:hybrid,)_with_states|
|-----|---|---|---|---|---|
|Akt pathway|2054|334|306|470|306|
|Bilirubin2_io|1680438|980307|980307| - |700275|
|Biohydrogenation_io|91|57|57|57|57|
|CD8 T cell differentiation|270|159|159|159|159|
|Chemical reaction network|81|78|78|78|78|
|Fujita|2054|334|306|470|306|
|Goodwin oscillator|138|138|138|514|514|
|HIV|694|694|342|1070|342|
|HIV2_io|439|439|231|245|231|
|JAK-STAT 1|553|475|475|501|475|
|LLW1987_io|705|705|705|273|1783|
|MAPK model (5 outputs bis)| - | - | - | - | - |
|MAPK model (5 outputs)| - | - | - | - | - |
|MAPK model (6 outputs)|595|595|595|595|595|
|Modified LV for testing|21|21|21|21|21|
|PK1|1886|1842|117|337|117|
|PK2| - | - | - | - | - |
|Pharm| - | - | - | - | - |
|QWWC|285|285| - |67|67|
|QY| - | - | - | - | - |
|SEAIJRC Covid model|1829|1829|1829|1829|1829|
|SEIR 34|187|187|187|103|103|
|SEIR 36 ref|210|210|210|210|210|
|SEIR_1_io|24752|23869|17585|17689|5927|
|SIR 19|55|55|55|55|55|
|SIR 21|55|55|55|55|55|
|SIR 24|29|29|95|29|46|
|SIR 6|57|21|21|21|21|
|SIRS forced|88|88|62|276|62|
|SIWR original|66|66|66|66|66|
|SIWR with extra output|66|66|66|66|66|
|SLIQR|790078|526676|226294|4217|5353|
|St|13645893|13645893|13645893| - | - |
|Treatment_io|139|88|88|88|78|
|TumorHu2019| - | - | - | - | - |
|TumorPillis2007| - | - | - | - | - |

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
