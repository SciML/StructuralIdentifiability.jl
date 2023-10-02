## Benchmark results

2023-09-13T19:02:49.963

  - Benchmarked function: `find_identifiable_functions`
  - Workers: 8
  - Timeout: 600 s

**All timings in seconds.**

For all models, `with_states = true`.

| Model                         | (:hybrid, 1) / Runtime | (:hybrid, 1) / Polynomial? | (:hybrid, 3) / Runtime | (:hybrid, 3) / Polynomial? |
|:----------------------------- |:---------------------- |:-------------------------- |:---------------------- |:-------------------------- |
| Akt pathway                   | 5.58                   | yes                        | 20.68                  | yes                        |
| Bilirubin2_io                 | 9.00                   | yes                        | 34.90                  | yes                        |
| Biohydrogenation_io           | 0.99                   | yes                        | 3.29                   | yes                        |
| Bruno2016                     | 0.61                   | yes                        | 2.02                   | yes                        |
| CD8 T cell differentiation    | 1.18                   | yes                        | 4.63                   | yes                        |
| CGV1990                       | 5.05                   | no                         | 10.20                  | no                         |
| Chemical reaction network     | 0.90                   | yes                        | 1.82                   | yes                        |
| Crauste_SI                    | 1.92                   | yes                        | 4.55                   | yes                        |
| Fujita                        | 9.32                   | yes                        | 21.76                  | yes                        |
| Goodwin oscillator            | 1.67                   | no                         | 3.83                   | no                         |
| HIV                           | 1.78                   | yes                        | 4.13                   | yes                        |
| HIV2_io                       | 14.63                  | no                         | 42.44                  | no                         |
| HighDimNonLin                 | 28.85                  | yes                        | 46.58                  | yes                        |
| KD1999                        | 2.21                   | no                         | 4.74                   | no                         |
| LLW1987_io                    | 0.88                   | no                         | 2.29                   | no                         |
| MAPK model (5 outputs)        | 58.75                  | yes                        | 66.04                  | yes                        |
| MAPK model (6 outputs)        | 13.01                  | yes                        | 21.02                  | yes                        |
| Modified LV for testing       | 0.33                   | yes                        | 0.66                   | yes                        |
| PK1                           | 2.62                   | yes                        | 7.38                   | yes                        |
| PK2                           | 168.70                 | yes                        | 168.62                 | yes                        |
| Pharm                         | 163.55                 | yes                        | 165.66                 | yes                        |
| Phosphorylation               | 0.98                   | yes                        | 1.97                   | yes                        |
| Pivastatin                    | 8.36                   | yes                        | 9.70                   | yes                        |
| QY                            | 73.06                  | no                         | 131.69                 | no                         |
| Ruminal lipolysis             | 0.47                   | yes                        | 1.18                   | yes                        |
| SEAIJRC Covid model           | 90.46                  | no                         | 95.61                  | no                         |
| SEIR 34                       | 1.57                   | yes                        | 4.15                   | yes                        |
| SEIR 36 ref                   | 2.69                   | yes                        | 4.83                   | yes                        |
| SEIR2T                        | 0.48                   | yes                        | 1.01                   | yes                        |
| SEIRT                         | 1.03                   | yes                        | 1.61                   | yes                        |
| SEIR_1_io                     | 1.41                   | no                         | 2.98                   | no                         |
| SEUIR                         | 1.13                   | no                         | 2.50                   | no                         |
| SIR 19                        | 0.70                   | yes                        | 1.73                   | yes                        |
| SIR 21                        | 0.69                   | yes                        | 2.02                   | yes                        |
| SIR 24                        | 0.73                   | yes                        | 1.60                   | yes                        |
| SIR 6                         | 0.49                   | yes                        | 1.29                   | yes                        |
| SIRS forced                   | 10.98                  | yes                        | 12.80                  | yes                        |
| SIWR original                 | 14.71                  | yes                        | 18.55                  | yes                        |
| SIWR with extra output        | 1.15                   | yes                        | 2.26                   | yes                        |
| SLIQR                         | 2.89                   | no                         | 8.08                   | no                         |
| St                            | 128.09                 | no                         | 600.0*                 | no                         |
| Transfection_4State           | 0.66                   | yes                        | 1.46                   | yes                        |
| Treatment_io                  | 2.17                   | yes                        | 4.76                   | yes                        |
| cLV1 (2o)                     | 2.09                   | yes                        | 4.11                   | yes                        |
| cholera                       | 1.14                   | yes                        | 1.91                   | yes                        |
| generalizedLoktaVolterra (1o) | 0.48                   | yes                        | 0.94                   | yes                        |
| generalizedLoktaVolterra (2o) | 0.44                   | yes                        | 0.93                   | yes                        |
| p53                           | 3.88                   | yes                        | 9.25                   | yes                        |

\* approx.

*Benchmarking environment:*

  - Total RAM (GiB): 15
  - Processor: Intel Xeon Processor (Icelake)
  - Julia version: 1.9.3

Versions of the dependencies:

  - Primes : 0.5.4
  - BenchmarkTools : 1.3.2
  - IterTools : 1.8.0
  - PrecompileTools : 1.2.0
  - Symbolics : 5.5.3
  - Combinatorics : 1.0.2
  - SymbolicUtils : 1.3.0
  - DataStructures : 0.18.15
  - Groebner : 0.4.3
  - ParamPunPam : 0.0.3
  - ModelingToolkit : 8.68.0
  - AbstractAlgebra : 0.31.1
  - MacroTools : 0.5.11
  - Nemo : 0.35.3
  - SpecialFunctions : 2.3.1
