@testset "Check identifiability of `ODESystem` object" begin
    @parameters a01 a21 a12
    @variables t x0(t) x1(t) y1(t) [output = true]
    D = Differential(t)

    eqs = [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1, y1 ~ x0]
    de = ODESystem(eqs, t, name = :Test)

    correct =
        Dict(a01 => :nonidentifiable, a21 => :nonidentifiable, a12 => :nonidentifiable)

    @test isequal(correct, assess_identifiability(de; measured_quantities = [y1 ~ x0]))
    @test isequal(correct, assess_identifiability(de; measured_quantities = [x0]))
    @test isequal(
        correct,
        assess_identifiability(de; measured_quantities = [(y1 ~ x0).rhs]),
    )

    # check identifiabile functions
    correct = [a01 * a12, a01 + a12 + a21]
    result = find_identifiable_functions(de, measured_quantities = [y1 ~ x0])
    @test isequal(Set(correct), Set(result))

    # --------------------------------------------------------------------------

    # check identifiabile functions
    @parameters V_m k_m k01 c
    @variables t x(t) y1(t) [output = true]
    D = Differential(t)

    eqs = [D(x) ~ (-V_m * x) / (k_m + x) + k01 * x, y1 ~ c * x]
    de = ODESystem(eqs, t, name = :Test)

    correct = [k01, c * k_m, V_m / k_m]
    result = find_identifiable_functions(de)
    @test isequal(Set(correct), Set(result))

    correct = [k01, c * x, k_m / x, V_m / x]
    result = find_identifiable_functions(de, with_states = true)
    @test isequal(Set(correct), Set(result))

    # --------------------------------------------------------------------------
    @parameters a01 a21 a12
    @variables t x0(t) x1(t) y1(t) [output = true]
    D = Differential(t)

    eqs = [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1, y1 ~ x0]
    de = ODESystem(eqs, t, name = :Test)

    correct =
        Dict(a01 => :nonidentifiable, a21 => :nonidentifiable, a12 => :nonidentifiable)

    @test isequal(correct, assess_identifiability(de))

    # --------------------------------------------------------------------------

    @parameters a01 a21 a12
    @variables t x0(t) x1(t) y1(t) [output = true]
    D = Differential(t)

    eqs = [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1, y1 ~ x0]
    de = ODESystem(eqs, t, name = :Test)
    funcs_to_check = [a01, a21, a12, a01 * a12, a01 + a12 + a21]
    correct = Dict(
        a12 => :nonidentifiable,
        a01 + a12 + a21 => :globally,
        a01 * a12 => :globally,
        a21 => :nonidentifiable,
        a01 => :nonidentifiable,
    )
    @test isequal(correct, assess_identifiability(de; funcs_to_check = funcs_to_check))

    # --------------------------------------------------------------------------

    @parameters a01 a21 a12
    @variables t x0(t) x1(t) y1(t)
    D = Differential(t)

    eqs = [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1]
    measured_quantities = [y1 ~ x0]
    de = ODESystem(eqs, t, name = :Test)
    funcs_to_check = [a01, a21, a12, a01 * a12, a01 + a12 + a21]
    correct = Dict(
        a12 => :nonidentifiable,
        a01 + a12 + a21 => :globally,
        a01 * a12 => :globally,
        a21 => :nonidentifiable,
        a01 => :nonidentifiable,
    )
    @test isequal(
        correct,
        assess_identifiability(
            de;
            measured_quantities = measured_quantities,
            funcs_to_check = funcs_to_check,
        ),
    )

    # --------------------------------------------------------------------------
    @parameters μ bi bw a χ γ k
    @variables t S(t) I(t) W(t) R(t) y(t)

    eqs = [
        D(S) ~ μ - bi * S * I - bw * S * W - μ * S + a * R,
        D(I) ~ bw * S * W + bi * S * I - (γ + μ) * I,
        D(W) ~ χ * (I - W),
        D(R) ~ γ * I - (μ + a) * R,
    ]
    de = ODESystem(eqs, t, name = :TestSIWR)
    measured_quantities = [y ~ k * I]
    # check all parameters (default)
    @test isequal(
        true,
        all(
            values(
                assess_local_identifiability(de; measured_quantities = measured_quantities),
            ),
        ),
    )

    # check specific parameters
    funcs_to_check = [μ, bi, bw, a, χ, γ, γ + μ, k, S, I, W, R]
    correct = Dict(f => true for f in funcs_to_check)
    @test isequal(
        correct,
        assess_local_identifiability(
            de;
            measured_quantities = measured_quantities,
            funcs_to_check = funcs_to_check,
        ),
    )

    # checking ME identifiability
    funcs_to_check = [μ, bi, bw, a, χ, γ, γ + μ, k]
    correct = Dict(f => true for f in funcs_to_check)
    @test isequal(
        (correct, 1),
        assess_local_identifiability(
            de;
            measured_quantities = measured_quantities,
            funcs_to_check = funcs_to_check,
            p = 0.99,
            type = :ME,
        ),
    )

    # checking identifiabile functions
    correct = [a, bw, χ, bi, k, γ, μ]
    result = find_identifiable_functions(de, measured_quantities = measured_quantities)
    @test isequal(Set(correct), Set(result))

    # --------------------------------------------------------------------------
    @parameters mu bi bw a xi gm k
    @variables t S(t) I(t) W(t) R(t) y(t) [output = true]

    eqs = [
        D(S) ~ mu - bi * S * I - bw * S * W - mu * S + a * R,
        D(I) ~ bw * S * W + bi * S * I - (gm + mu) * I,
        D(W) ~ xi * (I - W),
        D(R) ~ gm * I - (mu + a) * R,
        y ~ k * I,
    ]
    de = ODESystem(eqs, t, name = :TestSIWR)
    # check all parameters (default)
    @test isequal(true, all(values(assess_local_identifiability(de))))

    @test isequal(
        true,
        all(values(assess_local_identifiability(de; measured_quantities = [y ~ k * I]))),
    )

    # check specific parameters
    funcs_to_check = [mu, bi, bw, a, xi, gm, gm + mu, k, S, I, W, R]
    correct = Dict(f => true for f in funcs_to_check)
    @test isequal(
        correct,
        assess_local_identifiability(de; funcs_to_check = funcs_to_check),
    )

    # checking ME identifiability
    funcs_to_check = [mu, bi, bw, a, xi, gm, gm + mu, k]
    correct = Dict(f => true for f in funcs_to_check)
    @test isequal(
        (correct, 1),
        assess_local_identifiability(
            de;
            funcs_to_check = funcs_to_check,
            p = 0.99,
            type = :ME,
        ),
    )

    # --------------------------------------------------------------------------
    @parameters mu bi bw a xi gm k
    @variables t S(t) I(t) W(t) R(t) y(t)

    eqs = [
        D(S) ~ 2.0 * mu - bi * S * I - bw * S * W - mu * S + a * R,
        D(I) ~ bw * S * W + bi * S * I - (gm + mu) * I,
        D(W) ~ xi * (I - 0.6 * W),
        D(R) ~ gm * I - (mu + a) * R,
    ]
    de = ODESystem(eqs, t, name = :TestSIWR)
    measured_quantities = [y ~ 1.57 * I * k]
    funcs_to_check = [mu, bi, bw, a, xi, gm, mu, gm + mu, k, S, I, W, R]
    correct = Dict(f => true for f in funcs_to_check)
    @test isequal(
        correct,
        assess_local_identifiability(
            de;
            measured_quantities = measured_quantities,
            funcs_to_check = funcs_to_check,
        ),
    )

    # checking ME identifiability
    funcs_to_check = [bi, bw, a, xi, gm, mu, gm + mu, k]
    correct = Dict(f => true for f in funcs_to_check)
    @test isequal(
        (correct, 1),
        assess_local_identifiability(
            de;
            measured_quantities = measured_quantities,
            funcs_to_check = funcs_to_check,
            p = 0.99,
            type = :ME,
        ),
    )

    # ----------

    @parameters a01 a21 a12
    @variables t x0(t) x1(t) y1(t)
    D = Differential(t)
    using SpecialFunctions

    eqs = [
        D(x0) ~ -(a01 + a21) * SpecialFunctions.erfc(x0) + a12 * x1,
        D(x1) ~ a21 * x0 - a12 * x1,
    ]

    de = ODESystem(eqs, t, name = :Test)
    measured_quantities = [y1 ~ x0]
    funcs_to_check = [a01, a21, a12, a01 * a12, a01 + a12 + a21]
    correct = Dict(
        a01 => :nonidentifiable,
        a21 => :nonidentifiable,
        a12 => :nonidentifiable,
        a01 * a12 => :globally,
        a01 + a12 + a21 => :globally,
    )
    @test_throws ArgumentError assess_identifiability(
        de;
        measured_quantities = measured_quantities,
        funcs_to_check = funcs_to_check,
    )
    # ----------
    @parameters a b c
    @variables t x1(t) x2(t) y(t)
    D = Differential(t)

    eqs = [D(x1) ~ -a * x1 + x2 * b / (x1 + b / (c^2 - x2)), D(x2) ~ x2 * c^2 + x1]
    de = ODESystem(eqs, t, name = :Test)
    measured_quantities = [y ~ x2]
    correct = Dict(a => :globally, b => :globally, c => :locally)
    to_check = [a, b, c]
    @test isequal(
        correct,
        assess_identifiability(
            de;
            measured_quantities = measured_quantities,
            funcs_to_check = to_check,
        ),
    )

    # check identifiabile functions
    result = find_identifiable_functions(de, measured_quantities = measured_quantities)
    correct = [b, a, c^2]
    @test isequal(Set(result), Set(correct))

    # ----------
    @parameters a b
    @variables t c(t) x1(t) x2(t) y1(t) y2(t)
    D = Differential(t)

    eqs =
        [D(x1) ~ -a * x1 + x2 * b / (x1 + b / (c^2 - x2)), D(x2) ~ x2 * c^2 + x1, D(c) ~ 0]
    de = ODESystem(eqs, t, name = :Test)
    measured_quantities = [y1 ~ x2, y2 ~ c]
    correct = Dict(a => :globally, b => :globally)
    to_check = [a, b]
    @test isequal(
        correct,
        assess_identifiability(
            de;
            measured_quantities = measured_quantities,
            funcs_to_check = to_check,
        ),
    )

    #----------------------------------
    # Composable models test (from https://github.com/SciML/StructuralIdentifiability.jl/issues/162)
    @variables t
    function rabbits_creator(; name)
        ps = @parameters α = 1.5
        vars = @variables x(t) = 1.0 z(t) = 0.0 [input = true]
        D = Differential(t)
        equs = [D(x) ~ α^2 * x + z]

        ODESystem(equs, t, vars, ps; name = name)
    end

    function wolves_creator(; name)
        ps = @parameters δ = 3.0
        vars = @variables y(t) = 1.0 q(t) = 0.0 [input = true]
        D = Differential(t)
        equs = [D(y) ~ -δ * y + q]

        ODESystem(equs, t, vars, ps; name = name)
    end

    function lotka_volterra_creator(; name)
        @named wolves = wolves_creator()
        @named rabbits = rabbits_creator()

        ps = @parameters β = 1.0 γ = 1.0
        D = Differential(t)

        eqs = [rabbits.z ~ -β * wolves.y * rabbits.x, wolves.q ~ γ * wolves.y * rabbits.x]

        ModelingToolkit.compose(ODESystem(eqs, t, [], ps; name = name), wolves, rabbits)
    end

    @named ltk_mtk = lotka_volterra_creator()
    simp_ltk_mtk = structural_simplify(ltk_mtk)
    @unpack wolves₊δ, rabbits₊α, β, γ, wolves₊y, rabbits₊x = simp_ltk_mtk
    @unpack wolves₊y, rabbits₊x = simp_ltk_mtk
    @variables y(t)
    measured_quantities = [y ~ wolves₊y]
    result = assess_identifiability(simp_ltk_mtk, measured_quantities = measured_quantities)
    correct = Dict(
        rabbits₊α => :locally,
        γ => :nonidentifiable,
        β => :globally,
        wolves₊δ => :globally,
    )
    @test result == correct

    #----------------------------------

    @variables t, x(t), y(t), z(t), w(t)
    @parameters a
    @named sys = ODESystem([D(x) ~ a * y], t, [x], [a]; observed = [y ~ z, z ~ x])
    measured_quantities = [w ~ x]
    result = assess_identifiability(sys, measured_quantities = measured_quantities)
    @test result[a] == :globally

    result = find_identifiable_functions(sys, measured_quantities = measured_quantities)
    @test isequal(result, [a])

    #----------------------------------

    # Tensor definition case as reported in
    # https://github.com/SciML/StructuralIdentifiability.jl/issues/178
    @variables t, x(t)[1:2], y(t)[1:2]
    @parameters k1, k2

    eqs = [D(x[1]) ~ -k1 * x[2], D(x[2]) ~ -k2 * x[1]]

    sys = ODESystem(eqs, t, name = :example_vector)
    correct = Dict(k1 => true, k2 => true, x[1] => true, x[2] => true)
    @test assess_local_identifiability(sys, measured_quantities = [x[1], x[2]]) == correct
end
