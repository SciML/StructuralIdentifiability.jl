# Benchmark results

Benchmark models from [StructuralIdentifiability.jl/examples](https://github.com/SciML/StructuralIdentifiability.jl/tree/master/examples) and [https://github.com/Xabo-RB/Local-Global-Models](https://github.com/Xabo-RB/Local-Global-Models).

StructuralIdentifiability.jl on the `interface-id-funcs` branch.

All timings are in seconds.

## Aggregated results

For each model we report the runtime breakdown within StructuralIdentifiability.jl and the internal runtimes for Groebner bases in ParamPunPam.jl.

### StructuralIdentifiability.jl 

|Model|IO|global id.|extract funcs.|gen. ideal|inclusion|**ParamPunPam.jl**|total|
|-----|---|---|---|---|---|---|---|
SIWR with extra output|
|SIRS forced|
|Goodwin oscillator|
|Chemical reaction network|
|Akt pathway|
|LLW1987_io|
|HIV2_io|
|Biohydrogenation_io|
|Treatment_io|
|SLIQR|
|Fujita|
|NFkB| - | - | - | - | - | - | - |
|MAPK model (5 outputs)| 39.24 | 19.36 | - | 210.30 | - | - | - |
|Pharm| 14.33 | 25.71 | - | 109.98 | - | - | - |
|QWWC| 163.84 | 177.77 | - | - | - | - | - |
|SIWR original| 2.61 | 4.31 | - | 31.82 | - | - | - |
|MAPK model (6 outputs)| 6.36 | 10.97 | - | 39.02 | - | - | - |
|MAPK model (5 out. bis)| 412.44 | 631.71 | - | 5182.15 | - | - | - |
|SEAIJRC Covid model| 25.09 | 22.95 | - | 58.28 | - | - | - |
|Bilirubin2_io*| 0.01 | 0.05 | - | 0.67 | - | - | - |
|QY| 0.12 | 0.39 | - | 1.06 | - | - | - |
|St| 4.81 | 3.31 | - | 16.12 | - | - | - |
|SEIR_1_io*| 0.00 | 0.01 | - | 0.22 | - | - | - |
|HIV| 3.25 | 0.08 | - | 0.75 | - | - | - |
|CD8 T cell diff| 3.17 | 2.83 | - | 0.79 | - | - | - |

**- potentially an unlucky collision of evaluation points* 
