|Model|loc_time|glob_time|ioeq_time|wrnsk_time|rank_time|check_time|total|
|-----|---|---|---|---|---|---|---|
|Goodwin oscillator|0.03|0.11|0.03|0.06|0.00|0.02|0.15|
|Akt pathway|1.40|3.32|0.41|2.01|0.00|0.91|4.72|
|SIWR original|0.14|19.44|3.36|10.20|0.15|5.73|19.58|
|SEAIJRC Covid model|0.49|129.95|28.28|76.65|0.70|24.32|130.45|
|MAPK model (6 outputs)|5.39|25.46|1.99|11.89|0.00|11.57|30.85|
|Pharm|0.14|415.69|20.09|345.35|7.35|42.90|415.83|
|SIRS forced|0.04|27.38|0.80|25.42|0.23|0.93|27.42|
|SIWR with extra output|0.13|0.59|0.18|0.26|0.00|0.15|0.72|
|CD8 T cell differentiation|0.57|0.44|0.23|0.12|0.00|0.09|1.01|
|HIV|0.42|0.15|0.04|0.07|0.00|0.04|0.58|
|Chemical reaction network|0.26|0.22|0.07|0.11|0.00|0.04|0.48|
|Modified LV for testing|0.90|2.24|1.22|0.53|0.01|0.39|7.17|
|MAPK model (5 outputs)|5.00|55.00|27.41|13.86|0.00|13.72|60.00|

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
