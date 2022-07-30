|Model|loc_time|glob_time|ioeq_time|wrnsk_time|rank_time|check_time|total|
|-----|---|---|---|---|---|---|---|
|Goodwin oscillator|0.03|0.09|0.02|0.04|0.00|0.02|0.12|
|Akt pathway|0.64|4.37|0.18|0.95|0.00|3.23|5.01|
|SIWR original|0.08|13.08|2.97|2.86|0.14|7.10|13.16|
|SEAIJRC Covid model|0.07|76.34|27.58|16.23|0.63|31.90|76.41|
|MAPK model (6 outputs)|4.96|31.59|4.24|10.70|0.00|16.65|36.55|
|Pharm|0.05|100.86|19.30|39.82|5.64|36.10|100.90|
|SIRS forced|0.04|8.97|0.89|5.71|0.22|2.14|9.01|
|SIWR with extra output|0.09|0.67|0.44|0.09|0.00|0.13|0.76|
|CD8 T cell differentiation|0.48|0.20|0.02|0.09|0.00|0.09|0.69|
|HIV|0.11|0.41|0.04|0.06|0.00|0.31|0.52|
|MAPK model (5 outputs bis)|10.33|1665.97|461.87|214.31|0.05|989.75|1676.31|
|Chemical reaction network|0.04|0.14|0.05|0.05|0.00|0.04|0.18|
|Modified LV for testing|1.17|2.46|1.50|0.63|0.01|0.24|11.59|
|MAPK model (5 outputs)|3.91|98.06|59.62|9.97|0.00|28.46|101.97|

*Benchmarking environment:*

* Total RAM (Mb): 16384.0
* Processor: Intel(R) Core(TM) i5-8210Y CPU @ 1.60GHz
* Julia version: 1.7.2

Versions of the dependencies:

* Primes : 0.5.2
* BenchmarkTools : 1.3.1
* Singular : 0.10.1
* IterTools : 1.4.0
* DataStructures : 0.18.11
* Groebner : 0.2.5
* ModelingToolkit : 8.9.0
* AbstractAlgebra : 0.25.3
* MacroTools : 0.5.9
* Nemo : 0.30.0
* Hecke : 0.13.4
* SpecialFunctions : 2.1.4
* TestSetExtensions : 2.0.0
