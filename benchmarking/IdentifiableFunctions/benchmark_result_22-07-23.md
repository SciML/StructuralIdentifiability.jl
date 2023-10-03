## Benchmark results

2023-07-22T07:22:02.002

In this experiment:

  - Saturation variable had been moved to the last index (the smallest).
  - Globally identifiable parameters are not adjoined to the initial generators

| Model                      | io     | global id. | inclusion | inclusion Zp | ParamPunPam.jl | total  |
|:-------------------------- |:------ |:---------- |:--------- |:------------ |:-------------- |:------ |
| Akt pathway                | 2.67   | 0.00       | 1.40      | 0.89         | 3.82           | 9.38   |
| Bilirubin2_io              | 0.01   | 0.00       | 0.02      | 0.01         | 0.68           | 1.40   |
| Biohydrogenation_io        | 2.62   | 0.00       | 1.43      | 0.95         | 3.27           | 8.85   |
| CD8 T cell differentiation | 2.83   | 0.00       | 0.00      | 0.00         | 0.84           | 4.25   |
| Chemical reaction network  | 2.59   | 0.00       | 1.39      | 0.93         | 3.21           | 8.71   |
| Fujita                     | 2.65   | 0.00       | 1.40      | 0.89         | 3.83           | 9.38   |
| Goodwin oscillator         | 0.01   | 0.00       | 0.03      | 0.00         | 0.31           | 0.95   |
| HIV                        | 2.77   | 0.00       | 0.01      | 0.04         | 0.77           | 4.26   |
| HIV2_io                    | 2.39   | 0.00       | 0.00      | 0.00         | 0.81           | 3.83   |
| JAK-STAT 1                 | 10.29  | 0.00       | 1.65      | 1.17         | 4.34           | 18.13  |
| LLW1987_io                 | 0.00   | 0.00       | 1.38      | 0.91         | 3.26           | 6.12   |
| MAPK model (5 outputs bis) | -      | -          | -         | -            | -              | -      |
| MAPK model (5 outputs)     | 30.54  | 0.00       | 2.78      | 2.02         | 9.12           | 45.30  |
| MAPK model (6 outputs)     | 5.39   | 0.00       | 1.76      | 1.17         | 4.50           | 13.52  |
| Modified LV for testing    | 0.00   | 0.00       | 1.36      | 0.88         | 3.09           | 5.86   |
| NFkB                       | -      | -          | -         | -            | -              | -      |
| PK1                        | 2.39   | 0.00       | 0.01      | 0.00         | 0.74           | 3.74   |
| PK2                        | 14.79  | 0.00       | 18.27     | 4.53         | 10.63          | 51.04  |
| Pharm                      | 15.09  | 0.00       | 18.15     | 4.61         | 10.57          | 51.24  |
| QWWC                       | 162.29 | 0.00       | 119.75    | 26.37        | 39.70          | 368.93 |
| QY                         | 0.11   | 0.00       | 4.08      | 0.10         | 2.22           | 7.11   |
| SEAIJRC Covid model        | 22.86  | 0.00       | 15.00     | 2.99         | 11.24          | 54.77  |
| SEIR 34                    | 2.62   | 0.00       | 1.38      | 0.92         | 3.17           | 8.70   |
| SEIR 36 ref                | 2.57   | 0.00       | 0.01      | 0.01         | 0.70           | 3.88   |
| SEIR_1_io                  | 0.01   | 0.00       | 1.34      | 0.89         | 3.15           | 5.95   |
| SIR 19                     | 2.55   | 0.00       | 1.37      | 0.92         | 3.15           | 8.60   |
| SIR 21                     | 2.55   | 0.00       | 1.37      | 0.90         | 3.10           | 8.46   |
| SIR 24                     | 0.15   | 0.00       | 1.37      | 0.90         | 3.20           | 6.21   |
| SIR 6                      | 2.39   | 0.00       | 1.37      | 0.87         | 3.04           | 8.21   |
| SIRS forced                | 3.57   | 0.00       | 1.79      | 1.10         | 3.48           | 10.71  |
| SIWR original              | 2.50   | 0.00       | 3.18      | 0.71         | 1.85           | 9.20   |
| SIWR with extra output     | 3.25   | 0.00       | 0.09      | 0.02         | 0.40           | 4.38   |
| SLIQR                      | 0.01   | 0.00       | 1.41      | 0.92         | 3.61           | 6.56   |
| St                         | 4.21   | 0.00       | 3.47      | 0.70         | 15.73          | 26.19  |
| Treatment_io               | 2.39   | 0.00       | 1.38      | 0.90         | 3.11           | 8.33   |
| TumorHu2019                | -      | -          | -         | -            | -              | -      |
| TumorPillis2007            | -      | -          | -         | -            | -              | -      |

*Benchmarking environment:*

  - Total RAM (GiB): 2015
  - Processor: AMD EPYC 7702 64-Core Processor
  - Julia version: 1.9.1

Versions of the dependencies:

  - Primes : 0.5.3
  - BenchmarkTools : 1.3.2
  - IterTools : 1.8.0
  - PrecompileTools : 1.1.2
  - Symbolics : 5.5.0
  - SymbolicUtils : 1.0.5
  - DataStructures : 0.18.14
  - Groebner : 0.3.6
  - ParamPunPam : 0.0.1
  - ModelingToolkit : 8.60.0
  - AbstractAlgebra : 0.27.10
  - MacroTools : 0.5.10
  - Nemo : 0.32.7
  - SpecialFunctions : 2.3.0
