| Model                      | loc_time | glob_time | ioeq_time | wrnsk_time | rank_time | check_time | total   |
| -------------------------- | -------- | --------- | --------- | ---------- | --------- | ---------- | ------- |
| Goodwin oscillator         | 0.04     | 0.13      | 0.03      | 0.08       | 0.00      | 0.02       | 0.17    |
| Akt pathway                | 1.51     | 3.49      | 0.20      | 2.29       | 0.00      | 0.99       | 5.00    |
| SIWR original              | 0.12     | 17.89     | 3.05      | 9.35       | 0.15      | 5.33       | 18.01   |
| SEAIJRC Covid model        | 0.08     | 131.24    | 28.61     | 76.67      | 0.76      | 25.20      | 131.33  |
| MAPK model (6 outputs)     | 14.21    | 25.30     | 2.52      | 11.38      | 0.00      | 11.40      | 39.51   |
| Pharm                      | 0.04     | 405.70    | 18.42     | 337.13     | 7.22      | 42.94      | 405.74  |
| SIRS forced                | 0.04     | 30.21     | 0.97      | 27.98      | 0.24      | 1.01       | 30.25   |
| SIWR with extra output     | 0.13     | 0.54      | 0.16      | 0.23       | 0.00      | 0.14       | 0.67    |
| CD8 T cell differentiation | 0.36     | 0.44      | 0.02      | 0.11       | 0.00      | 0.31       | 0.80    |
| HIV                        | 0.09     | 0.49      | 0.37      | 0.08       | 0.00      | 0.04       | 0.57    |
| MAPK model (5 outputs bis) | 32.81    | 1051.67   | 397.40    | 225.49     | 0.05      | 428.72     | 1084.48 |
| Chemical reaction network  | 0.05     | 0.46      | 0.07      | 0.34       | 0.00      | 0.05       | 0.52    |
| Modified LV for testing    | 0.82     | 2.05      | 1.12      | 0.51       | 0.01      | 0.33       | 6.21    |
| MAPK model (5 outputs)     | 4.91     | 53.24     | 27.59     | 14.73      | 0.00      | 10.92      | 58.16   |

*Benchmarking environment:*

* Total RAM (Mb): 16384.0
* Processor: Intel(R) Core(TM) i5-8210Y CPU @ 1.60GHz
* Julia version: 1.6.1

Versions of the dependencies:

* Primes : 0.5.0
* Singular : 0.4.6
* IterTools : 1.3.0
* GroebnerBasis : 0.3.2
* AbstractAlgebra : 0.13.6
* MacroTools : 0.5.6
* Nemo : 0.20.1
* TestSetExtensions : 2.0.0
