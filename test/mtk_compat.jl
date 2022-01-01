@testset "Check identifiability of `ODESystem` object" begin
    @parameters a01 a21 a12 
    @variables t x0(t) x1(t) y1(t) [output=true]
    D = Differential(t)

    eqs = [
        D(x0) ~ -(a01 + a21) * x0 + a12 * x1,
        D(x1) ~ a21 * x0 - a12 * x1,
        y1 ~ x0
    ]
    de = ODESystem(eqs, t, name=:Test)

    correct = Dict(a01 => :nonidentifiable, a21 => :nonidentifiable, a12 => :nonidentifiable)

    @test isequal(correct, assess_identifiability(de, [y1~x0]))
    # --------------------------------------------------------------------------
    @parameters a01 a21 a12 
    @variables t x0(t) x1(t) y1(t) [output=true]
    D = Differential(t)

    eqs = [
        D(x0) ~ -(a01 + a21) * x0 + a12 * x1,
        D(x1) ~ a21 * x0 - a12 * x1,
        y1 ~ x0
    ]
    de = ODESystem(eqs, t, name=:Test)

    correct = Dict(a01 => :nonidentifiable, a21 => :nonidentifiable, a12 => :nonidentifiable)

    @test isequal(correct, assess_identifiability(de))

    # --------------------------------------------------------------------------

    @parameters a01 a21 a12 
    @variables t x0(t) x1(t) y1(t) [output=true]
    D = Differential(t)

    eqs = [
        D(x0) ~ -(a01 + a21) * x0 + a12 * x1,
        D(x1) ~ a21 * x0 - a12 * x1,
        y1 ~ x0
    ]
    de = ODESystem(eqs, t, name=:Test)
    funcs_to_check = [a01, a21, a12, a01 * a12, a01 + a12 + a21]
    correct = Dict(a12 => :nonidentifiable, a01 + a12 + a21 => :globally, a01 * a12 => :globally, a21 => :nonidentifiable, a01 => :nonidentifiable)
    @test isequal(correct, assess_identifiability(de, funcs_to_check))

    # --------------------------------------------------------------------------

    @parameters a01 a21 a12 
    @variables t x0(t) x1(t) y1(t)
    D = Differential(t)

    eqs = [
        D(x0) ~ -(a01 + a21) * x0 + a12 * x1,
        D(x1) ~ a21 * x0 - a12 * x1
    ]
    measured_quantities = [y1 ~ x0]
    de = ODESystem(eqs, t, name=:Test)
    funcs_to_check = [a01, a21, a12, a01 * a12, a01 + a12 + a21]
    correct = Dict(a12 => :nonidentifiable, a01 + a12 + a21 => :globally, a01 * a12 => :globally, a21 => :nonidentifiable, a01 => :nonidentifiable)
    @test isequal(correct, assess_identifiability(de, measured_quantities, funcs_to_check))

    # --------------------------------------------------------------------------
    @parameters mu bi bw a xi gm k
    @variables t S(t) I(t) W(t) R(t) y(t)

    eqs = [
        D(S) ~ mu - bi * S * I - bw * S * W - mu * S + a * R,
        D(I) ~ bw * S * W + bi * S * I - (gm + mu) * I,
        D(W) ~ xi * (I - W),
        D(R) ~ gm * I - (mu + a) * R
    ]
    de = ODESystem(eqs, t, name=:TestSIWR)
    measured_quantities = [y ~ k * I]
    # check all parameters (default)
    @test isequal(true, all(assess_local_identifiability(de, measured_quantities)))

    # check specific parameters
    funcs_to_check = [mu, bi, bw, a, xi, gm, gm + mu, k, S, I, W, R]
    correct = [true for _ in funcs_to_check]
    @test isequal(correct, assess_local_identifiability(de, measured_quantities, funcs_to_check))

    # checking ME identifiability
    funcs_to_check = [mu, bi, bw, a, xi, gm, gm + mu, k]
    correct = [true for _ in funcs_to_check] 
    @test isequal((correct, 1), assess_local_identifiability(de, measured_quantities, funcs_to_check, 0.99, :ME)) 

    # --------------------------------------------------------------------------
    @parameters mu bi bw a xi gm k
    @variables t S(t) I(t) W(t) R(t) y(t) [output=true]

    eqs = [
        D(S) ~ mu - bi * S * I - bw * S * W - mu * S + a * R,
        D(I) ~ bw * S * W + bi * S * I - (gm + mu) * I,
        D(W) ~ xi * (I - W),
        D(R) ~ gm * I - (mu + a) * R,
        y ~ k * I
    ]
    de = ODESystem(eqs, t, name=:TestSIWR)
    # check all parameters (default)
    @test isequal(true, all(assess_local_identifiability(de)))

    @test isequal(true, all(assess_local_identifiability(de, [y~k*I])))

    # check specific parameters
    funcs_to_check = [mu, bi, bw, a, xi, gm, gm + mu, k, S, I, W, R]
    correct = [true for _ in funcs_to_check]
    @test isequal(correct, assess_local_identifiability(de, funcs_to_check))

    # checking ME identifiability
    funcs_to_check = [mu, bi, bw, a, xi, gm, gm + mu, k]
    correct = [true for _ in funcs_to_check] 
    @test isequal((correct, 1), assess_local_identifiability(de, funcs_to_check, 0.99, :ME)) 
    
    # --------------------------------------------------------------------------
    @parameters mu bi bw a xi gm k
    @variables t S(t) I(t) W(t) R(t) y(t)

    eqs = [
        D(S) ~ 2.0 * mu - bi * S * I - bw * S * W - mu * S + a * R,
        D(I) ~ bw * S * W + bi * S * I - (gm + mu) * I,
        D(W) ~ xi * (I - 0.6 * W),
        D(R) ~ gm * I - (mu + a) * R,
    ]
    de = ODESystem(eqs, t, name=:TestSIWR)
    measured_quantities = [y ~ 1.57 * I * k]
    funcs_to_check = [mu, bi, bw, a, xi, gm, mu, gm + mu, k, S, I, W, R]
    correct = [true for _ in funcs_to_check]
    @test isequal(correct, assess_local_identifiability(de, measured_quantities, funcs_to_check))

    # checking ME identifiability
    funcs_to_check = [bi, bw, a, xi, gm, mu, gm + mu, k]
    correct = [true for _ in funcs_to_check] 
    @test isequal((correct, 1), assess_local_identifiability(de, measured_quantities, funcs_to_check, 0.99, :ME))
    
    # ----------

    @parameters a01 a21 a12 
    @variables t x0(t) x1(t) y1(t)
    D = Differential(t)
    using SpecialFunctions

    eqs = [
        D(x0) ~ -(a01 + a21) * SpecialFunctions.erfc(x0) + a12 * x1,
        D(x1) ~ a21 * x0 - a12 * x1
    ]

    de = ODESystem(eqs, t, name=:Test)
    measured_quantities = [y1 ~ x0]
    funcs_to_check = [a01, a21, a12, a01 * a12, a01 + a12 + a21]
    correct = Dict(a01 => :nonidentifiable, a21 => :nonidentifiable, a12 => :nonidentifiable, a01 * a12 => :globally, a01 + a12 + a21 => :globally)
    @test_throws ArgumentError assess_identifiability(de, measured_quantities,  funcs_to_check)
    # ----------
    @parameters a b c 
    @variables t x1(t) x2(t) y(t)
    D = Differential(t)

    eqs = [
        D(x1) ~ -a * x1 + x2 * b / (x1 + b / (c^2 - x2)),
        D(x2) ~ x2 * c^2 + x1
    ]
    de = ODESystem(eqs, t, name=:Test)
    measured_quantities = [y ~ x2]
    correct = Dict(a => :globally, b => :globally, c => :locally)
    to_check = [a, b, c]
    @test isequal(correct, assess_identifiability(de, measured_quantities, to_check))
end