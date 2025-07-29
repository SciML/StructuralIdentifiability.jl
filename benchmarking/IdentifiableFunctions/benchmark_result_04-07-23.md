# Benchmark results

Benchmark models from [StructuralIdentifiability.jl/examples](https://github.com/SciML/StructuralIdentifiability.jl/tree/master/examples) and [https://github.com/Xabo-RB/Local-Global-Models](https://github.com/Xabo-RB/Local-Global-Models).

StructuralIdentifiability.jl on the `interface-id-funcs` branch.

All timings are in seconds.

## Aggregated results

For each model we report the runtime breakdown within StructuralIdentifiability.jl and the internal runtimes for Groebner bases in ParamPunPam.jl.

### StructuralIdentifiability.jl

| Model | IO | global id. | extract funcs. | gen. ideal | inclusion | **ParamPunPam.jl** | total |
|:----- |:-- |:---------- |:-------------- |:---------- |:--------- |:------------------ |:----- |

SIWR with extra output|3.67|0.16|0.06|1.67|0.63|**372.89**|382.65|
|SIRS forced|4.36|0.85|0.16|4.21|1.92|**93.77**|107.31|
|Goodwin oscillator|0.02|0.11|0.0|0.78|0.06|**1.76**|2.77|
|Chemical reaction network|**3.11**|0.09|0.07|0.91|0.01|0.66|4.88|
|Akt pathway|3.26|3.91|0.54|1.48|0.61|**30.19**|43.45|
|LLW1987_io|0.0|0.01|0.0|0.0|0.03|**0.05**|0.12|
|HIV2_io|**2.98**|0.08|0.01|0.72|0.04|1.16|7.69|
|Biohydrogenation_io|**3.3**|0.06|0.01|0.7|0.01|0.29|4.41|
|Treatment_io|**3.06**|0.02|0.01|0.32|0.01|0.71|6.89|
|SLIQR|0.01|0.1|0.01|**0.92**|0.07|0.8|4.93|
|Fujita|3.49|4.24|0.63|1.58|0.63|**31.43**|45.49|
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

### Inside of ParamPunPam.jl

| Model | # points | discover shape | discover degrees | main loop | recover rationals |
|:----- |:-------- |:-------------- |:---------------- |:--------- |:----------------- |

SIWR with extra output|3906|2.07|6.27|364.09|0.0|
|SIRS forced|290|3.16|10.4|78.36|0.0|
|Goodwin oscillator|194|1.05|0.14|0.56|0.0|
|Chemical reaction network|194|0.11|0.1|0.44|0.0|
|Akt pathway|834|1.2|2.28|26.49|0.0|
|LLW1987_io|62|0.0|0.04|0.01|0.0|
|HIV2_io|194|0.84|0.08|0.24|0.0|
|Biohydrogenation_io|106|0.07|0.1|0.12|0.0|
|Treatment_io|610|0.1|0.07|0.53|0.0|
|SLIQR|514|0.07|0.07|0.66|0.0|
|Fujita|834|1.3|2.44|27.47|0.0|
|NFkB| - | - | - | - | - |
|MAPK model (5 outputs)| 608 | - | - | - | - |
|Pharm| 48 | - | - | - | - |
|QWWC| - | - | - | - | - |
|SIWR original| 2560 | - | - | - | - |
|MAPK model (6 outputs)| 3840 | - | - | - | - |
|MAPK model (5 out. bis)| 6 | - | - | - | - |
|SEAIJRC Covid model| 1600 | - | - | - | - |
|Bilirubin2_io| 229376 | - | - | - | - |
|QY| 360448 | - | - | - | - |
|St| 1344 | - | - | - | - |
|SEIR_1_io| 1179648 | - | - | - | - |
|HIV| 3407872 | - | - | - | - |
|CD8 T cell diff| 1310720 | - | - | - | - |

## Individual results

### SIWR with extra output

*Note: all parameters are globally identifiable.*

┌ Info: Given 141 functions in K(bi, gam, mu, bw, k, xi, a)[t, y1, y2, y3, y4, y5, y6, y7]

┌ Info: The shape of the basis is: 8 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y7, 1], [y6, 1], [y5, 1], [y4, 1], [y3, 1], [y2, 1], [y1, 1], [t, 1]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (0, 12)]]

┌ Info: Output summary:
│     Maximal interpolated degrees are: 1 for num. and 12 for den.
│     Maximal number of interpolated terms are: 1 for num. and 118 for den.

Interpolated functions:

```julia
AbstractAlgebra.Generic.Frac{Nemo.fmpq_mpoly}[
    a,
    xi,
    k,
    bw,
    mu,
    gam,
    bi,
    4 * bi^3 * mu^3 * k^3 * xi^3 - bi^2 * gam * mu^4 * k^3 * xi^2 - 6 * bi^2 * gam * mu^3 * k^3 * xi^3 + 2 * bi^2 * gam * mu^3 * k^3 * xi^2 * a - bi^2 * gam * mu^2 * k^3 * xi^4 - 6 * bi^2 * gam * mu^2 * k^3 * xi^3 * a + 4 * bi^2 * gam * mu^2 * k^3 * xi^2 * a^2 - bi^2 * mu^5 * k^3 * xi^2 - 6 * bi^2 * mu^4 * k^3 * xi^3 + 2 * bi^2 * mu^4 * k^3 * xi^2 * a + 12 * bi^2 * mu^3 * bw * k^3 * xi^3 - bi^2 * mu^3 * k^3 * xi^4 - 6 * bi^2 * mu^3 * k^3 * xi^3 * a + 4 * bi^2 * mu^3 * k^3 * xi^2 * a^2 + bi * gam^2 * mu^4 * k^3 * xi^2 - bi * gam^2 * mu^4 * k^3 * xi * a + 2 * bi * gam^2 * mu^3 * k^3 * xi^3 + bi * gam^2 * mu^3 * k^3 * xi^2 * a - 3 * bi * gam^2 * mu^3 * k^3 * xi * a^2 + bi * gam^2 * mu^2 * k^3 * xi^4 + 3 * bi * gam^2 * mu^2 * k^3 * xi^3 * a - bi * gam^2 * mu^2 * k^3 * xi^2 * a^2 - 3 * bi * gam^2 * mu^2 * k^3 * xi * a^3 + bi * gam^2 * mu * k^3 * xi^4 * a + bi * gam^2 * mu * k^3 * xi^3 * a^2 - bi * gam^2 * mu * k^3 * xi^2 * a^3 - bi * gam^2 * mu * k^3 * xi * a^4 + 2 * bi * gam * mu^5 * k^3 * xi^2 - 2 * bi * gam * mu^5 * k^3 * xi * a - 2 * bi * gam * mu^4 * bw * k^3 * xi^2 + 4 * bi * gam * mu^4 * k^3 * xi^3 + 2 * bi * gam * mu^4 * k^3 * xi^2 * a - 6 * bi * gam * mu^4 * k^3 * xi * a^2 - 12 * bi * gam * mu^3 * bw * k^3 * xi^3 + 4 * bi * gam * mu^3 * bw * k^3 * xi^2 * a + 2 * bi * gam * mu^3 * k^3 * xi^4 + 6 * bi * gam * mu^3 * k^3 * xi^3 * a - 2 * bi * gam * mu^3 * k^3 * xi^2 * a^2 - 6 * bi * gam * mu^3 * k^3 * xi * a^3 - 2 * bi * gam * mu^2 * bw * k^3 * xi^4 - 12 * bi * gam * mu^2 * bw * k^3 * xi^3 * a + 8 * bi * gam * mu^2 * bw * k^3 * xi^2 * a^2 + 2 * bi * gam * mu^2 * k^3 * xi^4 * a + 2 * bi * gam * mu^2 * k^3 * xi^3 * a^2 - 2 * bi * gam * mu^2 * k^3 * xi^2 * a^3 - 2 * bi * gam * mu^2 * k^3 * xi * a^4 + bi * mu^6 * k^3 * xi^2 - bi * mu^6 * k^3 * xi * a - 2 * bi * mu^5 * bw * k^3 * xi^2 + 2 * bi * mu^5 * k^3 * xi^3 + bi * mu^5 * k^3 * xi^2 * a - 3 * bi * mu^5 * k^3 * xi * a^2 - 12 * bi * mu^4 * bw * k^3 * xi^3 + 4 * bi * mu^4 * bw * k^3 * xi^2 * a + bi * mu^4 * k^3 * xi^4 + 3 * bi * mu^4 * k^3 * xi^3 * a - bi * mu^4 * k^3 * xi^2 * a^2 - 3 * bi * mu^4 * k^3 * xi * a^3 + 12 * bi * mu^3 * bw^2 * k^3 * xi^3 - 2 * bi * mu^3 * bw * k^3 * xi^4 - 12 * bi * mu^3 * bw * k^3 * xi^3 * a + 8 * bi * mu^3 * bw * k^3 * xi^2 * a^2 + bi * mu^3 * k^3 * xi^4 * a + bi * mu^3 * k^3 * xi^3 * a^2 - bi * mu^3 * k^3 * xi^2 * a^3 - bi * mu^3 * k^3 * xi * a^4 + gam^2 * mu^4 * bw * k^3 * xi^2 - gam^2 * mu^4 * bw * k^3 * xi * a + 2 * gam^2 * mu^3 * bw * k^3 * xi^3 + gam^2 * mu^3 * bw * k^3 * xi^2 * a - 3 * gam^2 * mu^3 * bw * k^3 * xi * a^2 + gam^2 * mu^2 * bw * k^3 * xi^4 + 3 * gam^2 * mu^2 * bw * k^3 * xi^3 * a - gam^2 * mu^2 * bw * k^3 * xi^2 * a^2 - 3 * gam^2 * mu^2 * bw * k^3 * xi * a^3 + gam^2 * mu * bw * k^3 * xi^4 * a + gam^2 * mu * bw * k^3 * xi^3 * a^2 - gam^2 * mu * bw * k^3 * xi^2 * a^3 - gam^2 * mu * bw * k^3 * xi * a^4 + 2 * gam * mu^5 * bw * k^3 * xi^2 - 2 * gam * mu^5 * bw * k^3 * xi * a - gam * mu^4 * bw^2 * k^3 * xi^2 + 4 * gam * mu^4 * bw * k^3 * xi^3 + 2 * gam * mu^4 * bw * k^3 * xi^2 * a - 6 * gam * mu^4 * bw * k^3 * xi * a^2 - 6 * gam * mu^3 * bw^2 * k^3 * xi^3 + 2 * gam * mu^3 * bw^2 * k^3 * xi^2 * a + 2 * gam * mu^3 * bw * k^3 * xi^4 + 6 * gam * mu^3 * bw * k^3 * xi^3 * a - 2 * gam * mu^3 * bw * k^3 * xi^2 * a^2 - 6 * gam * mu^3 * bw * k^3 * xi * a^3 - gam * mu^2 * bw^2 * k^3 * xi^4 - 6 * gam * mu^2 * bw^2 * k^3 * xi^3 * a + 4 * gam * mu^2 * bw^2 * k^3 * xi^2 * a^2 + 2 * gam * mu^2 * bw * k^3 * xi^4 * a + 2 * gam * mu^2 * bw * k^3 * xi^3 * a^2 - 2 * gam * mu^2 * bw * k^3 * xi^2 * a^3 - 2 * gam * mu^2 * bw * k^3 * xi * a^4 + mu^6 * bw * k^3 * xi^2 - mu^6 * bw * k^3 * xi * a - mu^5 * bw^2 * k^3 * xi^2 + 2 * mu^5 * bw * k^3 * xi^3 + mu^5 * bw * k^3 * xi^2 * a - 3 * mu^5 * bw * k^3 * xi * a^2 - 6 * mu^4 * bw^2 * k^3 * xi^3 + 2 * mu^4 * bw^2 * k^3 * xi^2 * a + mu^4 * bw * k^3 * xi^4 + 3 * mu^4 * bw * k^3 * xi^3 * a - mu^4 * bw * k^3 * xi^2 * a^2 - 3 * mu^4 * bw * k^3 * xi * a^3 + 4 * mu^3 * bw^3 * k^3 * xi^3 - mu^3 * bw^2 * k^3 * xi^4 - 6 * mu^3 * bw^2 * k^3 * xi^3 * a + 4 * mu^3 * bw^2 * k^3 * xi^2 * a^2 + mu^3 * bw * k^3 * xi^4 * a + mu^3 * bw * k^3 * xi^3 * a^2 - mu^3 * bw * k^3 * xi^2 * a^3 - mu^3 * bw * k^3 * xi * a^4,
]
```

## SIRS forced

*Note: the interpolated function that corresponds to parameters of degrees 4 can be expressed in terms of other functions of degrees <= 2.*

┌ Info: Given 2591 functions in K(nu, b1, b0, M, mu, g)[t, y1, y2, y3, y4, y5, y6]

┌ Info: The shape of the basis is: 6 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y6, 1], [y5, 1], [y3, 1], [y1, 1], [t, 1], [y4^2, 1]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (0, 4)], [(0, 0), (2, 0)]]

┌ Info: Output summary:
│     Maximal interpolated degrees are: 2 for num. and 4 for den.
│     Maximal number of interpolated terms are: 1 for num. and 15 for den.

Interpolated functions:

```julia
AbstractAlgebra.Generic.Frac{Nemo.fmpq_mpoly}[
    g,
    mu,
    b0,
    nu,
    M^2,
    nu^3 * mu + nu^3 * g + 2651 // 741 * nu^2 * mu^2 + 3242 // 741 * nu^2 * mu * g + 197 // 247 * nu^2 * g^2 + 46 // 247 * nu * M^2 * mu + 46 // 247 * nu * M^2 * g + 1161 // 247 * nu * mu^3 + 1421 // 247 * nu * mu^2 * g + 257 // 247 * nu * mu * g^2 - 3 // 247 * nu * g^3 + 46 // 247 * M^2 * mu^2 + 46 // 247 * M^2 * mu * g + 525 // 247 * mu^4 + 525 // 247 * mu^3 * g,
]
```

## Goodwin oscillator

┌ Info: Given 94 functions in K(b, alpha, c, gama, delta, sigma, beta)[t, y1, y2, y3, y4, y5, y6, y7]

┌ Info: The shape of the basis is: 6 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y6, 1], [y5, y7, 1], [y3, 1], [y1, 1], [t, 1], [y7^2, y7, 1]]

Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (1, 0)], [(0, 0), (0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (0, 6)], [(0, 0), (1, 0), (2, 0)]]

┌ Info: Output summary:
│     Maximal interpolated degrees are: 2 for num. and 6 for den.
│     Maximal number of interpolated terms are: 2 for num. and 5 for den.

Interpolated functions:

```julia
AbstractAlgebra.Generic.Frac{Nemo.fmpq_mpoly}[
    sigma,
    delta + beta,
    c,
    b,
    delta * beta,
    b^2 * c * sigma^2 + 1 // 3 * b^2 * c * sigma - 2 // 3 * b * c * delta * sigma^2 - 2 // 3 * b * c * sigma^2 * beta + c * delta * sigma^3 * beta,
]
```

## Chemical reaction network

*Note: all parameters are globally identifiable**

┌ Info: Given 78 functions in K(k5, k3, k4, k2, k6, k1)[t, y1, y2, y3, y4, y5, y6]

┌ Info: The shape of the basis is: 7 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y6, 1], [y5, 1], [y4, 1], [y3, 1], [y2, 1], [y1, 1], [t, 1]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (0, 7)]]

┌ Info: Output summary:
│     Maximal interpolated degrees are: 1 for num. and 7 for den.
│     Maximal number of interpolated terms are: 1 for num. and 6 for den.

Interpolated functions:

```julia
AbstractAlgebra.Generic.Frac{Nemo.fmpq_mpoly}[
    k1,
    k6,
    k2,
    k4,
    k3,
    k5,
    2 * k5^2 * k3^2 * k4 * k6 * k1 + 2 * k5^2 * k3 * k4 * k2 * k6 * k1 + 3 * k5 * k3^2 * k4^2 * k6 * k1 + 3 * k5 * k3 * k4^2 * k2 * k6 * k1 + k3^2 * k4^3 * k6 * k1 + k3 * k4^3 * k2 * k6 * k1,
]
```

## Akt pathway

**Too many variables: risk of overflow.**

**Result reported to be correct with probability 0.99**

┌ Info: Given 225 functions in K(a2, a3, reaction_1_k1, reaction_4_k1, reaction_9_k1, reaction_1_k2, reaction_7_k1, reaction_6_k1, a1, reaction_3_k1, reaction_5_k1, reaction_2_k1, EGFR_turnover, reaction_2_k2, reaction_5_k2, reaction_8_k1)[t, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16]

Info: The shape of the basis is: 13 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y16, 1], [y15, 1], [y14, 1], [y11, y12], [y10, 1], [y9, y12], [y8, 1], [y7, 1], [y4, 1], [y3, y5, y6, 1], [y2, y12], [y1, y12], [t*y12^8, 1]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 1)], [(0, 0), (1, 0)], [(0, 0), (1, 1)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (0, 0), (0, 0), (1, 0)], [(0, 0), (1, 1)], [(0, 0), (1, 1)], [(0, 0), (7, 15)]]

Warning: In Prime number approach the field order might be too small
│   ps = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59]
│   D = 15
│   order(K) = 4611686018427388039

Interpolated functions:

```julia
AbstractAlgebra.Generic.Frac{Nemo.fmpq_mpoly}[
    reaction_8_k1,
    reaction_5_k2,
    reaction_2_k2,
    reaction_5_k1 // reaction_2_k1,
    reaction_3_k1,
    a1 // reaction_2_k1,
    reaction_6_k1,
    reaction_7_k1,
    reaction_4_k1,
    reaction_1_k1 - reaction_9_k1 - reaction_1_k2,
    a3 // reaction_2_k1,
    a2 // reaction_2_k1,
    (1 // 16 * reaction_2_k1^7) // (a2^2 * a3^3 * reaction_1_k1 * reaction_7_k1 * reaction_6_k1^4 * a1 * reaction_3_k1 * reaction_5_k1 * reaction_8_k1 + a2^2 * a3^3 * reaction_1_k1 * reaction_7_k1 * reaction_6_k1^3 * a1 * reaction_3_k1 * reaction_5_k1 * reaction_8_k1^2 + a2^2 * a3^3 * reaction_1_k1 * reaction_6_k1^4 * a1 * reaction_3_k1^2 * reaction_5_k1 * reaction_8_k1 + a2^2 * a3^3 * reaction_1_k1 * reaction_6_k1^3 * a1 * reaction_3_k1^2 * reaction_5_k1 * reaction_8_k1^2 - a2^2 * a3^3 * reaction_4_k1 * reaction_7_k1 * reaction_6_k1^4 * a1 * reaction_3_k1 * reaction_5_k1 * reaction_8_k1 - a2^2 * a3^3 * reaction_4_k1 * reaction_7_k1 * reaction_6_k1^3 * a1 * reaction_3_k1 * reaction_5_k1 * reaction_8_k1^2 - a2^2 * a3^3 * reaction_4_k1 * reaction_6_k1^4 * a1 * reaction_3_k1^2 * reaction_5_k1 * reaction_8_k1 - a2^2 * a3^3 * reaction_4_k1 * reaction_6_k1^3 * a1 * reaction_3_k1^2 * reaction_5_k1 * reaction_8_k1^2 - a2^2 * a3^3 * reaction_9_k1 * reaction_7_k1 * reaction_6_k1^4 * a1 * reaction_3_k1 * reaction_5_k1 * reaction_8_k1 - a2^2 * a3^3 * reaction_9_k1 * reaction_7_k1 * reaction_6_k1^3 * a1 * reaction_3_k1 * reaction_5_k1 * reaction_8_k1^2 - a2^2 * a3^3 * reaction_9_k1 * reaction_6_k1^4 * a1 * reaction_3_k1^2 * reaction_5_k1 * reaction_8_k1 - a2^2 * a3^3 * reaction_9_k1 * reaction_6_k1^3 * a1 * reaction_3_k1^2 * reaction_5_k1 * reaction_8_k1^2 - a2^2 * a3^3 * reaction_1_k2 * reaction_7_k1 * reaction_6_k1^4 * a1 * reaction_3_k1 * reaction_5_k1 * reaction_8_k1 - a2^2 * a3^3 * reaction_1_k2 * reaction_7_k1 * reaction_6_k1^3 * a1 * reaction_3_k1 * reaction_5_k1 * reaction_8_k1^2 - a2^2 * a3^3 * reaction_1_k2 * reaction_6_k1^4 * a1 * reaction_3_k1^2 * reaction_5_k1 * reaction_8_k1 - a2^2 * a3^3 * reaction_1_k2 * reaction_6_k1^3 * a1 * reaction_3_k1^2 * reaction_5_k1 * reaction_8_k1^2),
]
```

## LLW1987_io

┌ Info: Given 13 functions in K(p2, p3, p4, p1)[t, y1, y2, y3, y4]

┌ Info: The shape of the basis is: 4 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y2, y4, 1], [t, 1], [y4^2, y4, 1], [y1*y3, 1]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (0, 0), (1, 0)], [(0, 0), (0, 3)], [(0, 0), (1, 0), (2, 0)], [(0, 0), (2, 0)]]

┌ Info: Output summary:
│     Maximal interpolated degrees are: 2 for num. and 3 for den.
│     Maximal number of interpolated terms are: 2 for num. and 2 for den.

Interpolated functions:

```julia
AbstractAlgebra.Generic.Frac{Nemo.fmpq_mpoly}[
    p3 + p1,
    p3 * p1,
    p2 * p4,
    p3^2 * p1 + p3 * p1^2,
]
```

## HIV2_io

┌ Info: Given 13 functions in K(b, c, q1, w1, k2, d, s, k1, w2, q2)[t, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10]

┌ Info: The shape of the basis is: 9 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y7, 1], [y6, 1], [y2, y4, y8, y9, 1], [y1, 1], [t, 1], [y5*y10, 1], [y3*y8, y4*y10, y8*y10, y10], [y4^2, y4*y8, y8^2, y4*y9, y8*y9, y9^2, y4, y8, y9, 1], [y9^3, y9^2, y9, 1]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (0, 0), (0, 0), (0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (0, 5)], [(0, 0), (2, 0)], [(0, 0), (0, 0), (0, 0), (2, 1)], [(0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (1, 0), (1, 0), (1, 0), (2, 0)], [(0, 0), (1, 0), (2, 0), (3, 0)]]

┌ Info: Output summary:
│     Maximal interpolated degrees are: 3 for num. and 5 for den.
│     Maximal number of interpolated terms are: 5 for num. and 3 for den.

Interpolated functions:

```julia
AbstractAlgebra.Generic.Frac{Nemo.fmpq_mpoly}[
    s,
    d,
    c + w1 + k1 + w2,
    b,
    k2 * q2,
    (q1 * k1 + w1 * q2 + k1 * q2) // q2,
    c * w1 + c * k1 + c * w2 + w1 * w2 + k1 * w2,
    c * w1 * w2 + c * k1 * w2,
    b * k2 * d * s * q2 - c * w1 * d * w2 - c * d * k1 * w2,
]
```

## Biohydrogenation_io

┌ Info: Given 41 functions in K(k5, k8, k9, k6, k10, k7)[t, y1, y2, y3, y4, y5, y6]

┌ Info: The shape of the basis is: 7 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y6, 1], [y4, 1], [y3, y5], [y2, y5, 1], [y1, 1], [t, y5], [y5^2, 1]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 1)], [(0, 0), (0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (0, 5)], [(0, 0), (2, 0)]]

┌ Info: Output summary:
│     Maximal interpolated degrees are: 2 for num. and 5 for den.
│     Maximal number of interpolated terms are: 2 for num. and 4 for den.
│

Interpolated functions:

```julia
AbstractAlgebra.Generic.Frac{Nemo.fmpq_mpoly}[
    k7,
    k6,
    k9 // k10,
    k8 + 1 // 2 * k10,
    k5,
    k10^2,
    k5 * k6 * k10^2 + 3 // 2 * k8^2 * k9 * k6 * k10 + 3 // 2 * k8 * k9 * k6 * k10^2 + 3 // 2 * k6 * k10^2 * k7,
]
```

## Treatment_io

┌ Info: Given 13 functions in K(b, nu, d, g, a)[t, y1, y2, y3, y4, y5]

┌ Info: The shape of the basis is: 8 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y2, y4, y5, 1], [y1, y4], [y3*y5, t, y3, y4, y5, 1], [t*y5, t, y3, 1], [y4^2, y4*y5, y5^2, y4, y5, 1], [y3*y4, y4, y5, 1], [t*y4, 1], [t^2, t*y3, y3^2]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (0, 0), (0, 0), (1, 0)], [(0, 0), (1, 1)], [(0, 0), (6, 0), (2, 0), (0, 0), (0, 0), (2, 0)], [(0, 0), (2, 0), (1, 3), (1, 3)], [(0, 0), (0, 0), (0, 0), (1, 0), (1, 0), (2, 0)], [(0, 0), (0, 0), (0, 0), (2, 0)], [(0, 0), (1, 3)], [(0, 0), (2, 6), (1, 9)]]

┌ Info: Output summary:
│     Maximal interpolated degrees are: 6 for num. and 9 for den.
│     Maximal number of interpolated terms are: 7 for num. and 10 for den.
│

Interpolated functions:

```julia
AbstractAlgebra.Generic.Frac{Nemo.fmpq_mpoly}[
    (1 // 4 * g) // (b * nu + b * d * g),
    (1 // 2 * g) // (b * nu + b * d * g),
    (1 // 4 * g) // (b^2 * nu^3 * d + 3 * b^2 * nu^2 * d^2 * g - b^2 * nu^2 * d * g - b^2 * nu^2 * d * a + 3 * b^2 * nu * d^3 * g^2 - 2 * b^2 * nu * d^2 * g^2 - 2 * b^2 * nu * d^2 * g * a + b^2 * d^4 * g^3 - b^2 * d^3 * g^3 - b^2 * d^3 * g^2 * a),
    nu + g + a,
    b // g,
    d * g - g - a,
    d * g - g - a,
    nu + d * g,
    (1 // 2 * nu + d * g - 1 // 2 * g - 1 // 2 * a) // (b * nu^2 * d + 2 * b * nu * d^2 * g - b * nu * d * g - b * nu * d * a + b * d^3 * g^2 - b * d^2 * g^2 - b * d^2 * g * a),
    nu * g + nu * a,
    4 * b * nu^2 * d + 8 * b * nu * d^2 * g - 4 * b * nu * d * g - 4 * b * nu * d * a + 4 * b * d^3 * g^2 - 4 * b * d^2 * g^2 - 4 * b * d^2 * g * a,
]
```

## SLIQR

*Note: it looks like more than 50% of runtime comes from functions that were not measured. Probably just moving the data around.*

┌ Info: Given 43 functions in K(b, e, Ninv, s, g, a)[t, y1, y2, y3, y4, y5, y6]

┌ Info: The shape of the basis is: 8 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y5, y6, 1], [y4, 1], [y3, 1], [y1, 1], [t, 1], [y6^2, y2, y6, 1], [y2*y6, y2, y6, 1], [y2^2, y2, y6, 1]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (0, 4)], [(0, 0), (4, 2), (3, 2), (5, 2)], [(0, 0), (3, 2), (2, 2), (3, 2)], [(0, 0), (7, 6), (3, 5), (7, 6)]]

┌ Info: Output summary:
│     Maximal interpolated degrees are: 7 for num. and 6 for den.
│     Maximal number of interpolated terms are: 14 for num. and 11 for den.

Interpolated functions:

```julia
AbstractAlgebra.Generic.Frac{Nemo.fmpq_mpoly}[
    g + a,
    s,
    Ninv,
    b,
    (e * a) // (e * s - s + a),
    (e * s * a) // (e * s - s + a),
    (e * s * a) // (e * s - s + a),
    (e * s^2 + e * s * g - s^2 - s * g + g * a + a^2) // (e * s - s + a),
    (e * s^2 * g - s^2 * g - s^2 * a + s * g * a + s * a^2) // (e * s - s + a),
    2 * e * Ninv * s * g + 2 * Ninv * s * a + 2 * Ninv * g * a,
    (e^2 * s * a - e^2 * a^2 - e * s * a + e * a^2) // (e^2 * s^2 * g - 2 * e * s^2 * g - e * s^2 * a + 2 * e * s * g * a + e * s * a^2 + s^2 * g + s^2 * a - 2 * s * g * a - 2 * s * a^2 + g * a^2 + a^3),
    (e^2 * s^2 * g - e * s^2 * g + 2 * e * s * g * a - s^2 * a - s * g * a + s * a^2 + g * a^2) // (e * s - s + a),
    (e^3 * s^2 * g * a - e^2 * s^3 * a - 2 * e^2 * s^2 * g * a + e^2 * s^2 * a^2 + 2 * e^2 * s * g * a^2 + e * s^3 * a + e * s^2 * g * a - e * s^2 * a^2 - 2 * e * s * g * a^2 + e * g * a^3) // (e^2 * s^3 * g - 2 * e * s^3 * g - e * s^3 * a + 2 * e * s^2 * g * a + e * s^2 * a^2 + s^3 * g + s^3 * a - 2 * s^2 * g * a - 2 * s^2 * a^2 + s * g * a^2 + s * a^3),
    (e^3 * s^3 * g - 2 * e^2 * s^3 * g - e^2 * s^3 * a + 3 * e^2 * s^2 * g * a + e^2 * s^2 * a^2 + e * s^3 * g - 4 * e * s^2 * g * a + 3 * e * s * g * a^2 + s^3 * a + s^2 * g * a - 2 * s^2 * a^2 - 2 * s * g * a^2 + s * a^3 + g * a^3) // (e^2 * s^3 * g - 2 * e * s^3 * g - e * s^3 * a + 2 * e * s^2 * g * a + e * s^2 * a^2 + s^3 * g + s^3 * a - 2 * s^2 * g * a - 2 * s^2 * a^2 + s * g * a^2 + s * a^3),
]
```

## Fujita

**Too many variables: risk of overflow.**

**Result reported to be correct with probability 0.99**

┌ Info: Given 225 functions in K(a2, a3, reaction_1_k1, reaction_4_k1, reaction_9_k1, reaction_1_k2, reaction_7_k1, reaction_6_k1, a1, reaction_3_k1, reaction_5_k1, reaction_2_k1, EGFR_turnover, reaction_2_k2, reaction_5_k2, reaction_8_k1)[t, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16]

┌ Info: The shape of the basis is: 13 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y16, 1], [y15, 1], [y14, 1], [y11, y12], [y10, 1], [y9, y12], [y8, 1], [y7, 1], [y4, 1], [y3, y5, y6, 1], [y2, y12], [y1, y12], [t*y12^8, 1]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 1)], [(0, 0), (1, 0)], [(0, 0), (1, 1)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (0, 0), (0, 0), (1, 0)], [(0, 0), (1, 1)], [(0, 0), (1, 1)], [(0, 0), (7, 15)]]

┌ Info: Output summary:
│     Maximal interpolated degrees are: 7 for num. and 15 for den.
│     Maximal number of interpolated terms are: 3 for num. and 16 for den.

┌ Warning: In Prime number approach the field order might be too small
│   ps = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59]
│   D = 15
│   order(K) = 4611686018427388039

Interpolated functions:

```julia
AbstractAlgebra.Generic.Frac{Nemo.fmpq_mpoly}[
    reaction_8_k1,
    reaction_5_k2,
    reaction_2_k2,
    reaction_5_k1 // reaction_2_k1,
    reaction_3_k1,
    a1 // reaction_2_k1,
    reaction_6_k1,
    reaction_7_k1,
    reaction_4_k1,
    reaction_1_k1 - reaction_9_k1 - reaction_1_k2,
    a3 // reaction_2_k1,
    a2 // reaction_2_k1,
    (1 // 16 * reaction_2_k1^7) // (a2^2 * a3^3 * reaction_1_k1 * reaction_7_k1 * reaction_6_k1^4 * a1 * reaction_3_k1 * reaction_5_k1 * reaction_8_k1 + a2^2 * a3^3 * reaction_1_k1 * reaction_7_k1 * reaction_6_k1^3 * a1 * reaction_3_k1 * reaction_5_k1 * reaction_8_k1^2 + a2^2 * a3^3 * reaction_1_k1 * reaction_6_k1^4 * a1 * reaction_3_k1^2 * reaction_5_k1 * reaction_8_k1 + a2^2 * a3^3 * reaction_1_k1 * reaction_6_k1^3 * a1 * reaction_3_k1^2 * reaction_5_k1 * reaction_8_k1^2 - a2^2 * a3^3 * reaction_4_k1 * reaction_7_k1 * reaction_6_k1^4 * a1 * reaction_3_k1 * reaction_5_k1 * reaction_8_k1 - a2^2 * a3^3 * reaction_4_k1 * reaction_7_k1 * reaction_6_k1^3 * a1 * reaction_3_k1 * reaction_5_k1 * reaction_8_k1^2 - a2^2 * a3^3 * reaction_4_k1 * reaction_6_k1^4 * a1 * reaction_3_k1^2 * reaction_5_k1 * reaction_8_k1 - a2^2 * a3^3 * reaction_4_k1 * reaction_6_k1^3 * a1 * reaction_3_k1^2 * reaction_5_k1 * reaction_8_k1^2 - a2^2 * a3^3 * reaction_9_k1 * reaction_7_k1 * reaction_6_k1^4 * a1 * reaction_3_k1 * reaction_5_k1 * reaction_8_k1 - a2^2 * a3^3 * reaction_9_k1 * reaction_7_k1 * reaction_6_k1^3 * a1 * reaction_3_k1 * reaction_5_k1 * reaction_8_k1^2 - a2^2 * a3^3 * reaction_9_k1 * reaction_6_k1^4 * a1 * reaction_3_k1^2 * reaction_5_k1 * reaction_8_k1 - a2^2 * a3^3 * reaction_9_k1 * reaction_6_k1^3 * a1 * reaction_3_k1^2 * reaction_5_k1 * reaction_8_k1^2 - a2^2 * a3^3 * reaction_1_k2 * reaction_7_k1 * reaction_6_k1^4 * a1 * reaction_3_k1 * reaction_5_k1 * reaction_8_k1 - a2^2 * a3^3 * reaction_1_k2 * reaction_7_k1 * reaction_6_k1^3 * a1 * reaction_3_k1 * reaction_5_k1 * reaction_8_k1^2 - a2^2 * a3^3 * reaction_1_k2 * reaction_6_k1^4 * a1 * reaction_3_k1^2 * reaction_5_k1 * reaction_8_k1 - a2^2 * a3^3 * reaction_1_k2 * reaction_6_k1^3 * a1 * reaction_3_k1^2 * reaction_5_k1 * reaction_8_k1^2),
]
```

## NFkB

**Process killed after 12 hours.**

Stuck on computing IO equations.

## MAPK model (5 outputs)

**Process killed after 12 hours.**

**Too many variables: risk of overflow.**

Used 608 interpolation points.

┌ Info: Given 527 functions in K(c0001, a10, gamma1000, alpha10, b00, beta11, c0111, beta10, alpha11, beta01, alpha01, gamma1100, c0011, c0010, b10, gamma1101, a00, b01, c1011, a01, gamma1110, gamma0100)[t, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16, y17, y18, y19, y20, y21, y22]

┌ Warning: In Prime number approach the field order might be too small
│   ps = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83]
│   D = 16
│   order(K) = 4611686018427388039

┌ Info: IO-equations computed in 39.24130316 seconds

┌ Info: Global identifiability assessed in 19.369218024 seconds

┌ Info: Differential ideal computed in 210.306699088 seconds

┌ Info: The shape of the basis is: 23 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y22, 1], [y21, 1], [y20, 1], [y19, 1], [y18, 1], [y17, 1], [y16, 1], [y15, 1], [y14, 1], [y13, 1], [y12, 1], [y11, 1], [y10, 1], [y9, 1], [y8, 1], [y7, 1], [y6, 1], [y5, 1], [y4, 1], [y3, 1], [y2, 1], [y1, 1], [t, 1]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (0, 16)]]

## Pharm

**Process killed after 12 hours.**

**Too many variables: risk of overflow.**

Used 48 interpolation points.

┌ Info: Given 3543 functions in K(a2, ka, n, b2, kc, b1, a1)[t, y1, y2, y3, y4, y5, y6, y7]

Warning: In Prime number approach the field order might be too small
│   ps = [2, 3, 5, 7, 11, 13, 17, 19]
│   D = 21
│   order(K) = 4611686018427388039

┌ Info: IO-equations computed in 14.333236304 seconds

Info: Global identifiability assessed in 25.719256061 seconds

┌ Info: Differential ideal computed in 109.984943125 seconds

┌ Info: The shape of the basis is: 8 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y7, 1], [y6, 1], [y5, 1], [y4, 1], [y3, 1], [y2, 1], [y1, 1], [t, 1]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (0, 21)]]

## QWWC

**Process killed.**

Process killed when producing the diff. ideal, maybe should look into this.

Info: IO-equations computed in 163.841701002 seconds

Info: Global identifiability assessed in 177.771413887 seconds

## SIWR original

**Process killed after 12 hours.**

**Too many variables: risk of overflow.**

Used 2560 interpolation points.

┌ Info: Given 778 functions in K(bi, gam, mu, bw, k, xi, a)[t, y1, y2, y3, y4, y5, y6, y7]

┌ Warning: In Prime number approach the field order might be too small
│   ps = [2, 3, 5, 7, 11, 13, 17, 19]
│   D = 17
│   order(K) = 4611686018427388039

┌ Info: IO-equations computed in 2.611795779 seconds

┌ Info: Global identifiability assessed in 4.31534692 seconds

┌ Info: Differential ideal computed in 31.826067139 seconds

┌ Info: The shape of the basis is: 8 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y7, 1], [y6, 1], [y5, 1], [y4, 1], [y3, 1], [y2, 1], [y1, 1], [t, 1]]

Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (0, 17)]]

## MAPK model (6 outputs)

**Process killed after 12 hours.**

**Too many variables: risk of overflow.**

Used 3840 interpolation points.

┌ Info: Given 117 functions in K(c0001, a10, gamma1000, alpha10, b00, beta11, c0111, beta10, alpha11, beta01, alpha01, gamma1100, c0011, c0010, b10, gamma1101, a00, b01, c1011, a01, gamma1110, gamma0100)[t, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16, y17, y18, y19, y20, y21, y22]

Warning: In Prime number approach the field order might be too small
│   ps = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83]
│   D = 12
│   order(K) = 4611686018427388039

┌ Info: IO-equations computed in 6.369373148 seconds

┌ Info: Global identifiability assessed in 10.970613039 seconds

┌ Info: Differential ideal computed in 39.025949794 seconds

┌ Info: The shape of the basis is: 23 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y22, 1], [y21, 1], [y20, 1], [y19, 1], [y18, 1], [y17, 1], [y16, 1], [y15, 1], [y14, 1], [y13, 1], [y12, 1], [y11, 1], [y10, 1], [y9, 1], [y8, 1], [y7, 1], [y6, 1], [y5, 1], [y4, 1], [y3, 1], [y2, 1], [y1, 1], [t, 1]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (0, 12)]]

## MAPK model (5 outputs bis)

**Process killed after 12 hours.**

**Too many variables: risk of overflow.**

Used 6 interpolation points.

Process killed before the parameter degrees are discovered.

┌ Info: IO-equations computed in 412.440120474 seconds

┌ Info: Global identifiability assessed in 631.717710978 seconds

┌ Info: Differential ideal computed in 5182.159844255 seconds

┌ Info: Given 1068 functions in K(c0001, a10, gamma1000, alpha10, b00, beta11, c0111, beta10, alpha11, beta01, alpha01, gamma1100, c0011, c0010, b10, gamma1101, a00, b01, c1011, a01, gamma1110, gamma0100)[t, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16, y17, y18, y19, y20, y21, y22]

┌ Info: The shape of the basis is: 23 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y22, 1], [y21, 1], [y20, 1], [y19, 1], [y18, 1], [y17, 1], [y16, 1], [y15, 1], [y14, 1], [y13, 1], [y12, 1], [y11, 1], [y10, 1], [y9, 1], [y8, 1], [y7, 1], [y6, 1], [y5, 1], [y4, 1], [y3, 1], [y2, 1], [y1, 1], [t, 1]]

## SEAIJRC Covid model

**Process killed after 12 hours.**

**Too many variables: risk of overflow.**

Used 1600 interpolation points.

┌ Info: Given 1470 functions in K(b, alpha, g2, k, g1, q, r)[t, y1, y2, y3, y4, y5, y6, y7]

┌ Warning: In Prime number approach the field order might be too small
│   ps = [2, 3, 5, 7, 11, 13, 17, 19]
│   D = 17
│   order(K) = 4611686018427388039

┌ Info: IO-equations computed in 25.096497945 seconds

┌ Info: Global identifiability assessed in 22.950870056 seconds

┌ Info: Differential ideal computed in 58.286382806 seconds

┌ Info: The shape of the basis is: 8 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y5, 1], [y4, 1], [y3, 1], [y2, 1], [y1, 1], [y6*y7, y6, y7], [t*y6^2, t*y7^2, t*y6, y6^2, t*y7, y6, 1], [t*y7^3, t*y7^2, t*y6, t*y7, y6, 1]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (0, 0), (2, 1)], [(0, 0), (4, 2), (2, 1), (4, 15), (4, 2), (5, 15), (6, 15)], [(0, 0), (0, 0), (1, 2), (0, 0), (5, 17), (4, 15)]]

## Bilirubin2_io

**Process killed after 12 hours.**

Potentially a bug:

Does not finish even with 229,376 points.
Just 5k points are enough for interpolating 7 variables of degree 3. Something probably went wrong.

┌ Info: Given 10 functions in K(k01, k31, k21, k12, k13, k14, k41)[t, y1, y2, y3, y4, y5, y6, y7]

┌ Info: IO-equations computed in 0.012185094 seconds

┌ Info: Global identifiability assessed in 0.053605773 seconds

┌ Info: Differential ideal computed in 0.675370293 seconds

┌ Info: The shape of the basis is: 13 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y4, y5, y6, 1], [y2, y3, y7, 1], [y1, 1], [t, 1], [y7^2, y6, y7, 1], [y6*y7, y6, y7, 1], [y3*y7, y3, y5, y7, 1], [y6^2, y6, y7, 1], [y5*y6, y3, y5, y6, 1], [y3*y6, y5*y7, y3, y5, y6, y7, 1], [y5^2, y3, y5, y7, 1], [y3*y5, y5*y7, y3, y5, y7, 1], [y3^2, y3, y5, y6, 1]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (0, 0), (0, 0), (1, 0)], [(0, 0), (0, 0), (0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (0, 3)], [(0, 0), (0, 0), (1, 0), (2, 0)], [(0, 0), (1, 0), (1, 0), (2, 0)], [(0, 0), (1, 0), (0, 0), (1, 0), (2, 0)], [(0, 0), (1, 0), (0, 0), (2, 0)], [(0, 0), (0, 0), (1, 0), (1, 0), (2, 0)], [(0, 0), (0, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0)], [(0, 0), (0, 0), (1, 0), (0, 0), (2, 0)], [(0, 0), (0, 0), (1, 0), (1, 0), (1, 0), (2, 0)], [(0, 0), (1, 0), (0, 0), (0, 0), (2, 0)]]

## QY

**Process killed after 12 hours.**

Used 360448 interpolation points.

┌ Info: Given 68 functions in K(siga1, beta_SI, phi, alpa, M, Mar, Ks, beta, siga2, beta_SA)[t, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10]

┌ Info: IO-equations computed in 0.128419756 seconds

┌ Info: Global identifiability assessed in 0.397637536 seconds

┌ Info: Differential ideal computed in 1.066933301 seconds

┌ Info: The shape of the basis is: 15 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y7, 1], [y6, 1], [y4, 1], [y3, y5, y9, 1], [y1, y5, y9, 1], [t, y5, y9, y10, 1], [y10^2, y5, y9, y10, 1], [y5*y10, y9*y10, y5, y9, y10, 1], [y2*y10, y9*y10, y2, y9, y10, 1], [y9^2, y2, y5, y9, y10, 1], [y5*y9, y10, 1], [y2*y9, y9*y10, y10, 1], [y5^2, y2, y5, 1], [y2*y5, y2, y5, 1], [y2^2, y2, y5, 1]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0), (1, 0), (2, 0)], [(0, 0), (0, 0), (0, 0), (1, 0)], [(0, 0), (4, 11), (4, 11), (4, 12), (5, 11)], [(0, 0), (7, 4), (7, 4), (5, 3), (8, 4)], [(0, 0), (0, 0), (3, 1), (3, 1), (1, 0), (4, 1)], [(0, 0), (4, 3), (3, 1), (5, 3), (5, 3), (6, 3)], [(0, 0), (3, 3), (4, 3), (1, 0), (3, 3), (2, 0)], [(0, 0), (3, 3), (4, 2)], [(0, 0), (0, 0), (4, 3), (4, 2)], [(0, 0), (3, 3), (4, 3), (4, 2)], [(0, 0), (4, 3), (5, 4), (4, 2)], [(0, 0), (7, 6), (6, 7), (7, 5)]]

## St

**Process killed after 12 hours.**

**Too many variables: risk of overflow.**

Used 1344 interpolation points.

┌ Info: Given 42 functions in K(e, rR, dr, d, g, r, a, T, Dd)[t, y1, y2, y3, y4, y5, y6, y7, y8, y9]

┌ Warning: In Prime number approach the field order might be too small
│   ps = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
│   D = 21
│   order(K) = 4611686018427388039

┌ Info: IO-equations computed in 4.813259583 seconds

┌ Info: Global identifiability assessed in 3.31958208 seconds

┌ Info: Differential ideal computed in 16.128728705 seconds

┌ Info: The shape of the basis is: 26 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y9, 1], [y8, 1], [y1, y2, y3, y4, y5, y6, y7, 1], [y3*y5, y4*y5, y3, y4, y7], [y2*y5, y5*y6, y2, y3, y4, y7, 1], [y4^2, y4*y6, y2*y7, y3*y7, y4*y7, y6*y7, y3, y4, y7], [y3*y4, y3*y6, y2*y7, y3*y7, y3, y4, y7], [y2*y4, y3*y6, y3, y4], [y3^2, y3*y6, y4*y6, y2*y7, y3*y7, y4*y7, y6*y7, y7^2, y3, y4, y7], [y2*y3, y3*y6, y4*y6, y2*y7, y3*y7, y4*y7, y6*y7, y7^2, y3, y4, y7], [y2^2, y3*y6, y4*y6, y2*y7, y3*y7, y4*y7, y6*y7, y7^2, y2, y3, y4, y7, 1], [y4*y5*y7, y5*y6*y7, y4*y5, y2*y7, y3*y7, y5*y7, y3, y4, y7], [y5*y6^2, y5*y6*y7, y5*y7^2, y2*y6, y3*y6, y4*y6, y5*y6, y2*y7, y3*y7, y4*y7, y5*y7, y6*y7, y7^2, y2, y3, y4, y5, y6, y7, 1], [y4*y5*y6, y5*y6*y7, y5*y7^2, y4*y5, y3*y6, y4*y6, y2*y7, y3*y7, y4*y7, y5*y7, y6*y7, y7^2, y3, y4, y7], [t*y7^5, y4*y5^3, y5^3*y6, y4*y5^2, y5^3, y5^2*y6, y4*y5, y5^2, y2*y6, y3*y6, y4*y6, y5*y6, y6^2, y2*y7, y3*y7, y4*y7, y6*y7, y2, y3, y4, y5, y6, y7, 1], [t*y4*y7^4, y4*y5^3, y5^3*y6, y4*y5^2, y5^3, y5^2*y6, y4*y5, y5^2, y2*y6, y3*y6, y4*y6, y5*y6, y6^2, y2*y7, y3*y7, y4*y7, y6*y7, y2, y3, y4, y5, y6, y7, 1], [t*y3*y7^4, y4*y5^3, y5^3*y6, y4*y5^2, y5^3, y5^2*y6, y4*y5, y5^2, y2*y6, y3*y6, y5*y6, y2*y7, y3*y7, y2, y3, y4, y5, y6, y7, 1], [t*y4*y6*y7^3, t*y2*y7^4, t*y6*y7^4, t*y3*y7^3, t*y4*y7^3, t*y7^4, y4*y5^3, y5^3*y6, y4*y5^2, y5^3, y5^2*y6, y4*y5, y5^2, y5*y6, y2, y3, y4, y5, y6, 1], [t*y3*y6*y7^3, t*y2*y7^4, t*y3*y7^3, t*y4*y7^3, t*y7^4, y4*y5^3, y5^3*y6, y4*y5^2, y5^3, y5^2*y6, y4*y5, y5^2, y5*y6, y2, y3, y5, 1], [t*y4*y6^2*y7^2, t*y2*y6*y7^3, t*y6^2*y7^3, t*y3*y6*y7^2, t*y4*y6*y7^2, t*y2*y7^3, t*y3*y7^3, t*y4*y7^3, t*y6*y7^3, y4*y5^3, y5^3*y6, t*y3*y7^2, t*y4*y7^2, t*y7^3, y4*y5^2, y5^3, y5^2*y6, y4*y5, y5^2, y5*y6, y2, y3, y4, y5, y6, 1], [t*y3*y6^2*y7^2, t*y2*y6*y7^3, t*y3*y6*y7^2, t*y4*y6*y7^2, t*y2*y7^3, t*y6*y7^3, y4*y5^3, y5^3*y6, t*y3*y7^2, t*y4*y7^2, t*y7^3, y4*y5^2, y5^3, y5^2*y6, y4*y5, y5^2, y5*y6, y2, y3, y5, 1], [t*y4*y6^3*y7, t*y2*y6^2*y7^2, t*y6^3*y7^2, t*y3*y6^2*y7, t*y4*y6^2*y7, t*y2*y6*y7^2, t*y3*y6*y7^2, t*y4*y6*y7^2, t*y6^2*y7^2, t*y2*y7^3, t*y3*y7^3, t*y4*y7^3, t*y6*y7^3, y4*y5^3, y5^3*y6, t*y3*y6*y7, t*y4*y6*y7, t*y2*y7^2, t*y3*y7^2, t*y4*y7^2, t*y6*y7^2, t*y7^3, y4*y5^2, y5^3, y5^2*y6, t*y3*y7, t*y4*y7, t*y7^2, y4*y5, y5^2, y5*y6, y2, y3, y4, y5, y6, 1], [t*y3*y6^3*y7, t*y2*y6^2*y7^2, t*y3*y6^2*y7, t*y4*y6^2*y7, t*y2*y6*y7^2, t*y6^2*y7^2, y4*y5^3, y5^3*y6, t*y3*y6*y7, t*y4*y6*y7, t*y2*y7^2, t*y3*y7^2, t*y4*y7^2, t*y6*y7^2, y4*y5^2, y5^3, y5^2*y6, t*y3*y7, t*y4*y7, t*y7^2, y4*y5, y5^2, y5*y6, y2, y3, y5, 1], [t*y4*y6^4, t*y2*y6^3*y7, t*y6^4*y7, t*y3*y6^3, t*y4*y6^3, t*y2*y6^2*y7, t*y3*y6^2*y7, t*y4*y6^2*y7, t*y6^3*y7, t*y2*y6*y7^2, t*y3*y6*y7^2, t*y4*y6*y7^2, t*y6^2*y7^2, t*y2*y7^3, t*y3*y7^3, t*y4*y7^3, t*y6*y7^3, y4*y5^3, y5^3*y6, t*y3*y6^2, t*y4*y6^2, t*y2*y6*y7, t*y3*y6*y7, t*y4*y6*y7, t*y6^2*y7, t*y2*y7^2, t*y3*y7^2, t*y4*y7^2, t*y6*y7^2, t*y7^3, y4*y5^2, y5^3, t*y3*y6, t*y4*y6, y5^2*y6, t*y2*y7, t*y3*y7, t*y4*y7, t*y6*y7, t*y7^2, t*y3, t*y4, y4*y5, y5^2, y5*y6, t*y7, y2, y3, y4, y5, y6, 1], [t*y3*y6^4, t*y2*y6^3*y7, t*y3*y6^3, t*y4*y6^3, t*y2*y6^2*y7, t*y6^3*y7, y4*y5^3, y5^3*y6, t*y3*y6^2, t*y4*y6^2, t*y2*y6*y7, t*y3*y6*y7, t*y4*y6*y7, t*y6^2*y7, t*y2*y7^2, t*y3*y7^2, t*y4*y7^2, t*y6*y7^2, y4*y5^2, y5^3, t*y3*y6, t*y4*y6, y5^2*y6, t*y2*y7, t*y3*y7, t*y4*y7, t*y6*y7, t*y7^2, t*y3, t*y4, y4*y5, y5^2, y5*y6, t*y7, y2, y3, y5, 1], [t*y5*y6*y7^4, t*y5*y6*y7^3, t*y2*y7^4, t*y5*y7^4, y4*y5^4, y5^4*y6, t*y5*y6*y7^2, t*y2*y7^3, t*y3*y7^3, t*y4*y7^3, t*y5*y7^3, t*y7^4, y4*y5^3, y5^4, y5^3*y6, t*y5*y6*y7, t*y2*y7^2, t*y3*y7^2, t*y4*y7^2, t*y5*y7^2, t*y7^3, t*y4*y5, y4*y5^2, y5^3, y5^2*y6, t*y2*y7, t*y3*y7, t*y4*y7, t*y5*y7, t*y7^2, t*y3, t*y4, y4*y5, y5^2, y5*y6, t*y7, y2, y3, y4, y5, y6, 1]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (0, 0), (1, 0), (1, 0), (0, 0), (0, 0), (1, 0), (2, 0)], [(0, 0), (0, 0), (2, 1), (2, 1), (2, 1)], [(0, 0), (0, 0), (1, 0), (3, 1), (3, 1), (3, 1), (3, 0)], [(0, 0), (0, 1), (0, 1), (0, 0), (0, 0), (0, 1), (3, 2), (3, 2), (3, 2)], [(0, 0), (0, 1), (0, 1), (0, 0), (3, 2), (3, 2), (3, 2)], [(0, 0), (0, 0), (1, 0), (1, 0)], [(0, 0), (2, 3), (0, 1), (2, 3), (0, 0), (0, 0), (0, 1), (0, 0), (5, 4), (3, 2), (5, 4)], [(0, 0), (2, 2), (0, 0), (2, 2), (1, 0), (1, 0), (0, 0), (1, 0), (3, 2), (1, 0), (4, 2)], [(0, 0), (3, 2), (1, 0), (3, 2), (2, 0), (2, 0), (1, 0), (2, 0), (2, 0), (6, 3), (4, 1), (6, 3), (4, 0)], [(0, 0), (0, 1), (3, 2), (1, 1), (1, 0), (3, 2), (5, 3), (5, 3), (5, 3)], [(0, 0), (1, 0), (4, 2), (1, 0), (5, 3), (5, 3), (2, 0), (4, 2), (6, 3), (6, 3), (2, 0), (5, 3), (6, 3), (3, 0), (7, 3), (7, 3), (4, 0), (3, 0), (7, 3), (5, 0)], [(0, 0), (0, 0), (3, 2), (3, 1), (4, 3), (4, 3), (3, 2), (5, 3), (5, 3), (3, 1), (4, 3), (5, 3), (7, 4), (7, 4), (7, 4)], [(0, 0), (11, 19), (10, 19), (12, 19), (12, 19), (11, 19), (13, 19), (13, 19), (0, 11), (1, 11), (1, 11), (12, 19), (0, 11), (1, 11), (2, 11), (2, 11), (1, 11), (14, 20), (16, 21), (4, 12), (14, 19), (2, 11), (4, 12), (16, 20)], [(0, 0), (11, 19), (10, 19), (12, 19), (12, 19), (11, 19), (13, 19), (13, 19), (0, 11), (1, 11), (1, 11), (12, 19), (0, 11), (1, 11), (2, 11), (2, 11), (1, 11), (14, 20), (16, 21), (4, 12), (14, 19), (2, 11), (4, 12), (16, 20)], [(0, 0), (11, 19), (10, 19), (12, 19), (12, 19), (11, 19), (13, 19), (13, 19), (0, 11), (1, 11), (12, 19), (1, 11), (2, 11), (14, 20), (16, 21), (3, 11), (14, 19), (1, 11), (4, 12), (16, 20)], [(0, 0), (0, 0), (0, 0), (3, 1), (3, 1), (3, 1), (12, 19), (11, 19), (13, 19), (13, 19), (12, 19), (14, 19), (14, 19), (13, 19), (14, 19), (15, 19), (3, 10), (15, 19), (2, 10), (16, 19)], [(0, 0), (0, 0), (3, 1), (3, 1), (3, 1), (12, 19), (11, 19), (13, 19), (13, 19), (12, 19), (14, 19), (14, 19), (13, 19), (14, 19), (15, 19), (15, 19), (16, 19)], [(0, 0), (0, 0), (0, 0), (3, 1), (3, 1), (2, 0), (3, 0), (3, 0), (3, 1), (13, 19), (12, 19), (5, 1), (5, 1), (5, 1), (14, 19), (14, 19), (13, 19), (15, 19), (15, 19), (14, 19), (15, 19), (16, 19), (4, 10), (16, 19), (3, 10), (17, 19)], [(0, 0), (0, 0), (3, 1), (3, 1), (1, 0), (3, 1), (13, 19), (12, 19), (4, 1), (5, 1), (4, 1), (14, 19), (14, 19), (13, 19), (15, 19), (15, 19), (14, 19), (15, 19), (16, 19), (16, 19), (17, 19)], [(0, 0), (0, 0), (0, 0), (3, 1), (3, 1), (2, 0), (3, 0), (3, 0), (3, 1), (3, 0), (4, 0), (4, 0), (3, 0), (14, 19), (13, 19), (5, 1), (5, 1), (4, 0), (5, 0), (5, 0), (5, 1), (5, 0), (15, 19), (15, 19), (14, 19), (7, 1), (7, 1), (7, 1), (16, 19), (16, 19), (15, 19), (16, 19), (17, 19), (5, 10), (17, 19), (4, 10), (18, 19)], [(0, 0), (0, 0), (3, 1), (3, 1), (1, 0), (3, 1), (14, 19), (13, 19), (4, 1), (5, 1), (4, 0), (5, 0), (5, 0), (5, 1), (15, 19), (15, 19), (14, 19), (7, 1), (7, 1), (7, 1), (16, 19), (16, 19), (15, 19), (16, 19), (17, 19), (17, 19), (18, 19)], [(0, 0), (0, 0), (0, 0), (3, 1), (3, 1), (2, 0), (3, 0), (3, 0), (3, 1), (3, 0), (4, 0), (4, 0), (3, 0), (4, 0), (5, 0), (5, 0), (4, 0), (15, 19), (14, 19), (5, 1), (5, 1), (4, 0), (5, 0), (5, 0), (5, 1), (5, 0), (6, 0), (6, 0), (5, 0), (6, 0), (16, 19), (16, 19), (7, 1), (7, 1), (15, 19), (6, 0), (8, 1), (8, 1), (7, 1), (8, 1), (9, 1), (9, 1), (17, 19), (17, 19), (16, 19), (9, 1), (17, 19), (18, 19), (6, 10), (18, 19), (5, 10), (19, 19)], [(0, 0), (0, 0), (3, 1), (3, 1), (1, 0), (3, 1), (15, 19), (14, 19), (4, 1), (5, 1), (4, 0), (5, 0), (5, 0), (5, 1), (5, 0), (6, 0), (6, 0), (5, 0), (16, 19), (16, 19), (7, 1), (7, 1), (15, 19), (6, 0), (8, 1), (8, 1), (7, 1), (8, 1), (9, 1), (9, 1), (17, 19), (17, 19), (16, 19), (9, 1), (17, 19), (18, 19), (18, 19), (19, 19)], [(0, 0), (3, 2), (1, 0), (3, 1), (12, 19), (11, 19), (6, 4), (4, 2), (5, 2), (5, 2), (6, 3), (5, 2), (13, 19), (13, 19), (12, 19), (9, 6), (7, 4), (8, 4), (8, 4), (9, 5), (8, 4), (12, 7), (14, 19), (14, 19), (13, 19), (10, 6), (11, 6), (11, 6), (12, 7), (11, 6), (14, 8), (14, 8), (15, 19), (15, 19), (14, 19), (14, 8), (15, 19), (16, 19), (4, 10), (16, 19), (3, 10), (17, 19)]]

## SEIR_1_io

**Process killed after 12 hours.**

Potentially a bug:

Does not finish even with 1,179,648 points.
Just 30k points are enough for interpolating 4 variables of degree 5. Something probably went wrong.

┌ Info: Given 19 functions in K(beta, gamma, v, psi)[t, y1, y2, y3, y4]

┌ Info: IO-equations computed in 0.009532897 seconds

┌ Info: Global identifiability assessed in 0.016378922 seconds

┌ Info: Differential ideal computed in 0.228718016 seconds

┌ Info: The shape of the basis is: 5 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y3, y4, 1], [y2, 1], [y1, y4], [t, y4, 1], [y4^2, y4, 1]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (1, 0), (2, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 1)], [(0, 0), (1, 5), (2, 5)], [(0, 0), (2, 1), (1, 1)]]

## HIV

**Process killed after 12 hours.**

Used 3407872 interpolation points.

┌ Info: Given 87 functions in K(b, c, h, lm, d, k, u, q, a, beta)[t, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10]

┌ Info: IO-equations computed in 3.2581641 seconds

┌ Info: Global identifiability assessed in 0.082135079 seconds

┌ Info: Differential ideal computed in 0.753892455 seconds

┌ Info: The shape of the basis is: 10 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y9, 1], [y7, 1], [y5, 1], [y4, y8], [y3, 1], [y1, 1], [t, 1], [y2*y8, y6*y10], [y6*y8*y10, 1], [y6^2*y10^2, y2]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 1)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (0, 7)], [(0, 0), (2, 2)], [(0, 0), (3, 0)], [(0, 0), (4, 1)]]

## CD8 T cell diff.

**Process killed after 12 hours.**

Used 1310720 interpolation points.

┌ Info: IO-equations computed in 3.172457507 seconds

┌ Info: Global identifiability assessed in 2.830197433 seconds

┌ Info: Differential ideal computed in 0.79689821 seconds

┌ Info: Given 41 functions in K(rho_P, mu_EE, delta_NE, mu_LE, delta_EL, mu_N, mu_M, delta_LM, mu_PE, mu_PL, mu_LL, mu_P, rho_E)[t, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13]

┌ Info: The shape of the basis is: 13 polynomials with monomials
│   state.shape = Vector{Nemo.gfp_mpoly}[[y12, 1], [y11, 1], [y10, 1], [y9, 1], [y8, 1], [y7, 1], [y6, 1], [y5, 1], [y4, 1], [y3, y13], [y2, 1], [y1, y13], [t*y13, 1]]

┌ Info: The total degrees in the coefficients
│   state.param_degrees = [[(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 0)], [(0, 0), (1, 1)], [(0, 0), (1, 0)], [(0, 0), (1, 1)], [(0, 0), (1, 7)]]

* * *

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
