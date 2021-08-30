@testset "Check identifiability of `ODESystem` object" begin
    @parameters a01 a21 a12 
    @variables t x0(t) x1(t) y1(t) [output = true]
    D = Differential(t)

    eqs = [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1, y1 ~ x0]

    de = ODESystem(eqs, t, name=:Test)
    funcs_to_check = [a01, a21, a12, a01 * a12, a01 + a12 + a21]
    correct = [:nonidentifiable, :nonidentifiable, :nonidentifiable, :globally, :globally]
    @test isequal(correct, assess_identifiability(de, funcs_to_check))

    # --------------------------------------------------------------------------
    @parameters μ bi bw a ξ γ k
    @variables t S(t) I(t) W(t) R(t) y(t) [output = true]

    eqs = [
        D(S) ~ μ - bi * S * I - bw * S * W - μ * S + a * R,
        D(I) ~ bw * S * W + bi * S * I - (γ + μ) * I,
        D(W) ~ ξ * (I - W),
        D(R) ~ γ * I - (μ + a) * R,
        y ~ k * I
    ]
    de = ODESystem(eqs, t, name=:TestSIWR)
    # check all parameters (default)
    @test isequal(true, all(assess_local_identifiability(de)))

    # check specific parameters
    funcs_to_check = [μ, bi, bw, a, ξ, γ, μ, γ + μ, k, S, I, W, R]
    correct = [true for _ in funcs_to_check]
    @test isequal(correct, assess_local_identifiability(de, funcs_to_check))

    # checking ME identifiability
    funcs_to_check = [μ, bi, bw, a, ξ, γ, μ, γ + μ, k]
    correct = [true for _ in funcs_to_check] 
    @test isequal((correct, 1), assess_local_identifiability(de, funcs_to_check, 0.99, :ME)) 

    # --------------------------------------------------------------------------
    @parameters μ bi bw a ξ γ k
    @variables t S(t) I(t) W(t) R(t) y(t) [output = true]

    eqs = [
        D(S) ~ 2.0 * μ - bi * S * I - bw * S * W - μ * S + a * R,
        D(I) ~ bw * S * W + bi * S * I - (γ + μ) * I,
        D(W) ~ ξ * (I - 0.6 * W),
        D(R) ~ γ * I - (μ + a) * R,
        y ~ 1.57 * I * k
    ]
    de = ODESystem(eqs, t, name=:TestSIWR)
    funcs_to_check = [μ, bi, bw, a, ξ, γ, μ, γ + μ, k, S, I, W, R]
    correct = [true for _ in funcs_to_check]
    @test isequal(correct, assess_local_identifiability(de, funcs_to_check))

    # checking ME identifiability
    funcs_to_check = [μ, bi, bw, a, ξ, γ, μ, γ + μ, k]
    correct = [true for _ in funcs_to_check] 
    @test isequal((correct, 1), assess_local_identifiability(de, funcs_to_check, 0.99, :ME))
    
    # ----------

    @parameters a01 a21 a12 
    @variables t x0(t) x1(t) y1(t) [output = true]
    D = Differential(t)
    using SpecialFunctions

    eqs = [
        D(x0) ~ -(a01 + a21) * SpecialFunctions.erfc(x0) + a12 * x1,
        D(x1) ~ a21 * x0 - a12 * x1,
        y1 ~ x0
    ]

    de = ODESystem(eqs, t, name=:Test)
    
    funcs_to_check = [a01, a21, a12, a01 * a12, a01 + a12 + a21]
    correct = [:nonidentifiable, :nonidentifiable, :nonidentifiable, :globally, :globally]
    @test_throws ArgumentError assess_identifiability(de,  funcs_to_check)
    # ----------
    @parameters a b c 
    @variables t x1(t) x2(t) y(t) [output = true]
    D = Differential(t)

    eqs = [D(x1) ~ -a * x1 + x2 * b / (x1 + b / (c^2 - x2)), D(x2) ~ x2 * c^2 + x1, y ~ x2]
    de = ODESystem(eqs, t, name=:Test)
    correct = [:globally, :globally, :locally]
    to_check = [a, b, c]
    @test isequal(correct, assess_identifiability(de, to_check))

end