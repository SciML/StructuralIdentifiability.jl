@testset "Check identifiability of `ODESystem` object" begin
    @parameters a01 a21 a12 
    @variables t x0(t) x1(t) y1(t)
    D = Differential(t)

    eqs = [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1]

    de = ODESystem(eqs, t, [x0, x1], [a01, a21, a12], name=:Test)
    output_eqs = [y1~x0]
    inputs = []
    funcs_to_check = [a01, a21, a12, a01 * a12, a01 + a12 + a21]
    correct = [:nonidentifiable, :nonidentifiable, :nonidentifiable, :globally, :globally]
    @test isequal(correct, assess_identifiability(de, output_eqs, inputs, funcs_to_check))

    #--------------------------------------------------------------------------
    @parameters μ bi bw a ξ γ k
    @variables t S(t) I(t) W(t) R(t) y(t)


    eqs = [
        D(S) ~ μ - bi * S * I - bw * S * W - μ * S + a * R,
        D(I) ~ bw * S * W + bi * S * I - (γ + μ) * I,
        D(W) ~ ξ * (I - W),
        D(R) ~ γ * I - (μ + a) * R
    ]
    de = ODESystem(eqs, t, [S, I, W, R], [μ, bi, bw, a, ξ, γ, k], name=:TestSIWR)
    output_eqs = [y ~ k * I]
    funcs_to_check = [μ, bi, bw, a, ξ, γ, μ, γ + μ, k, S, I, W, R]
    correct = [true for _ in funcs_to_check]
    inputs = []
    @test isequal(correct, assess_local_identifiability(de, output_eqs, inputs, funcs_to_check))

    # checking ME identifiability
    funcs_to_check = [μ, bi, bw, a, ξ, γ, μ, γ + μ, k]
    correct = [true for _ in funcs_to_check] 
    @test isequal(, assess_local_identifiability(de, output_eqs, inputs, funcs_to_check, 0.99, :ME)) 


    #--------------------------------------------------------------------------

end