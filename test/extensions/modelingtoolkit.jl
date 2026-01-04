if GROUP == "All" || GROUP == "ModelingToolkitSIExt"
    @testset "eval_at_nemo" begin
        using ModelingToolkitBase, Symbolics
        using Nemo

        @independent_variables t
        @parameters a01 a21 a12
        @variables x0 x1

        ring, (a, b, c, x, y) = QQ["a", "b", "c", "x", "y"]

        nemo = StructuralIdentifiability.eval_at_nemo(
            x0 + x1 * a01^2 + x1^20 * (a21 + a12),
            Dict(x0 => x, x1 => y, a01 => a, a21 => b, a12 => c),
        )
        @test nemo == x + y * a^2 + y^20 * (b + c)
    end

    @testset "Check identifiability of `System` object" begin
        using ModelingToolkitBase
        using ModelingToolkitBase: parameters
        using Symbolics

        @independent_variables t
        @parameters a01 a21 a12
        @variables x0(t) x1(t) y1(t) [output = true]
        D = Differential(t)

        eqs = [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1, y1 ~ x0]
        de = System(eqs, t, name = :Test)

        correct = OrderedDict(
            a01 => :nonidentifiable,
            a21 => :nonidentifiable,
            a12 => :nonidentifiable,
            x0 => :globally,
            x1 => :nonidentifiable,
        )

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
        @independent_variables t
        @parameters V_m k_m k01 c
        @variables x(t) y1(t) [output = true]
        D = Differential(t)

        eqs = [D(x) ~ (-V_m * x) / (k_m + x) + k01 * x, y1 ~ c * x]
        de = System(eqs, t, name = :Test)

        correct = [k01, c * k_m, V_m * c]
        result = find_identifiable_functions(de)
        @test isequal(Set(correct), Set(result))

        correct = [k01, c * x, k_m * c, V_m * c]
        result = find_identifiable_functions(de, with_states = true)
        @test isequal(Set(correct), Set(result))

        # --------------------------------------------------------------------------
        @independent_variables t
        @parameters a01 a21 a12
        @variables x0(t) x1(t) y1(t) [output = true]
        D = Differential(t)

        eqs = [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1, y1 ~ x0]
        de = System(eqs, t, name = :Test)

        correct = OrderedDict(
            a01 => :nonidentifiable,
            a21 => :nonidentifiable,
            a12 => :nonidentifiable,
            x0 => :globally,
            x1 => :nonidentifiable,
        )

        @test isequal(correct, assess_identifiability(de))

        # --------------------------------------------------------------------------

        @independent_variables t
        @parameters a01 a21 a12
        @variables x0(t) x1(t) y1(t) [output = true]
        D = Differential(t)

        eqs = [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1, y1 ~ x0]
        de = System(eqs, t, name = :Test)
        funcs_to_check = [a01, a21, a12, a01 * a12, a01 + a12 + a21]
        correct = OrderedDict(
            a01 => :nonidentifiable,
            a21 => :nonidentifiable,
            a12 => :nonidentifiable,
            a01 * a12 => :globally,
            a01 + a12 + a21 => :globally,
        )
        @test isequal(correct, assess_identifiability(de; funcs_to_check = funcs_to_check))

        # --------------------------------------------------------------------------

        @independent_variables t
        @parameters a01 a21 a12
        @variables x0(t) x1(t) y1(t)
        D = Differential(t)

        eqs = [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1]
        measured_quantities = [y1 ~ x0]
        de = System(eqs, t, name = :Test)
        funcs_to_check = [a01, a21, a12, a01 * a12, a01 + a12 + a21]
        correct = OrderedDict(
            a01 => :nonidentifiable,
            a21 => :nonidentifiable,
            a12 => :nonidentifiable,
            a01 * a12 => :globally,
            a01 + a12 + a21 => :globally,
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
        @independent_variables t
        @parameters μ bi bw a χ γ k
        @variables S(t) I(t) W(t) R(t) y(t)

        eqs = [
            D(S) ~ μ - bi * S * I - bw * S * W - μ * S + a * R,
            D(I) ~ bw * S * W + bi * S * I - (γ + μ) * I,
            D(W) ~ χ * (I - W),
            D(R) ~ γ * I - (μ + a) * R,
        ]
        de = System(eqs, t, name = :TestSIWR)
        measured_quantities = [y ~ k * I]
        # check all parameters (default)
        @test isequal(
            true,
            all(
                values(
                    assess_local_identifiability(
                        de;
                        measured_quantities = measured_quantities,
                    ),
                ),
            ),
        )

        # check specific parameters
        funcs_to_check = [μ, bi, bw, a, χ, γ, γ + μ, k, S, I, W, R]
        correct = OrderedDict(f => true for f in funcs_to_check)
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
        correct = OrderedDict(f => true for f in funcs_to_check)
        @test isequal(
            (correct, 1),
            assess_local_identifiability(
                de;
                measured_quantities = measured_quantities,
                funcs_to_check = funcs_to_check,
                prob_threshold = 0.99,
                type = :ME,
            ),
        )

        # checking identifiabile functions
        correct = [a, bw, χ, bi, k, γ, μ]
        result = find_identifiable_functions(de, measured_quantities = measured_quantities)
        @test isequal(Set(correct), Set(result))

        # --------------------------------------------------------------------------
        @independent_variables t
        @parameters mu bi bw a xi gm k
        @variables S(t) I(t) W(t) R(t) y(t) [output = true]

        eqs = [
            D(S) ~ mu - bi * S * I - bw * S * W - mu * S + a * R,
            D(I) ~ bw * S * W + bi * S * I - (gm + mu) * I,
            D(W) ~ xi * (I - W),
            D(R) ~ gm * I - (mu + a) * R,
            y ~ k * I,
        ]
        de = System(eqs, t, name = :TestSIWR)
        # check all parameters (default)
        @test isequal(true, all(values(assess_local_identifiability(de))))

        @test isequal(
            true,
            all(
                values(assess_local_identifiability(de; measured_quantities = [y ~ k * I])),
            ),
        )

        # check specific parameters
        funcs_to_check = [mu, bi, bw, a, xi, gm, gm + mu, k, S, I, W, R]
        correct = OrderedDict(f => true for f in funcs_to_check)
        @test isequal(
            correct,
            assess_local_identifiability(de; funcs_to_check = funcs_to_check),
        )

        # checking ME identifiability
        funcs_to_check = [mu, bi, bw, a, xi, gm, gm + mu, k]
        correct = OrderedDict(f => true for f in funcs_to_check)
        @test isequal(
            (correct, 1),
            assess_local_identifiability(
                de;
                funcs_to_check = funcs_to_check,
                prob_threshold = 0.99,
                type = :ME,
            ),
        )

        # --------------------------------------------------------------------------
        @independent_variables t
        @parameters mu bi bw a xi gm k
        @variables S(t) I(t) W(t) R(t) y(t)

        eqs = [
            D(S) ~ 2.0 * mu - bi * S * I - bw * S * W - mu * S + a * R,
            D(I) ~ bw * S * W + bi * S * I - (gm + mu) * I,
            D(W) ~ xi * (I - 0.6 * W),
            D(R) ~ gm * I - (mu + a) * R,
        ]
        de = System(eqs, t, name = :TestSIWR)
        measured_quantities = [y ~ 1.57 * I * k]
        funcs_to_check = [mu, bi, bw, a, xi, gm, mu, gm + mu, k, S, I, W, R]
        correct = OrderedDict(f => true for f in funcs_to_check)
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
        correct = OrderedDict(f => true for f in funcs_to_check)
        @test isequal(
            (correct, 1),
            assess_local_identifiability(
                de;
                measured_quantities = measured_quantities,
                funcs_to_check = funcs_to_check,
                prob_threshold = 0.99,
                type = :ME,
            ),
        )

        # ----------

        @independent_variables t
        @parameters a01 a21 a12
        @variables x0(t) x1(t) y1(t)
        D = Differential(t)
        using SpecialFunctions

        eqs = [
            D(x0) ~ -(a01 + a21) * SpecialFunctions.erfc(x0) + a12 * x1,
            D(x1) ~ a21 * x0 - a12 * x1,
        ]

        de = System(eqs, t, name = :Test)
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
        @independent_variables t
        @parameters a b c
        @variables x1(t) x2(t) y(t)
        D = Differential(t)

        eqs = [D(x1) ~ -a * x1 + x2 * b / (x1 + b / (c^2 - x2)), D(x2) ~ x2 * c^2 + x1]
        de = System(eqs, t, name = :Test)
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
        @independent_variables t
        @parameters a b
        @variables c(t) x1(t) x2(t) y1(t) y2(t)
        D = Differential(t)

        eqs = [
            D(x1) ~ -a * x1 + x2 * b / (x1 + b / (c^2 - x2)),
            D(x2) ~ x2 * c^2 + x1,
            D(c) ~ 0,
        ]
        de = System(eqs, t, name = :Test)
        measured_quantities = [y1 ~ x2, y2 ~ c]
        correct = OrderedDict(a => :globally, b => :globally)
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
        @independent_variables t
        function rabbits_creator(; name)
            ps = @parameters α = 1.5
            vars = @variables x(t) = 1.0 z(t) = 0.0 [input = true]
            D = Differential(t)
            equs = [D(x) ~ α^2 * x + z]

            System(equs, t, vars, ps; name = name)
        end

        function wolves_creator(; name)
            ps = @parameters δ = 3.0
            vars = @variables y(t) = 1.0 q(t) = 0.0 [input = true]
            D = Differential(t)
            equs = [D(y) ~ -δ * y + q]

            System(equs, t, vars, ps; name = name)
        end

        function lotka_volterra_creator(; name)
            @named wolves = wolves_creator()
            @named rabbits = rabbits_creator()

            ps = @parameters β = 1.0 γ = 1.0
            D = Differential(t)

            eqs =
                [rabbits.z ~ -β * wolves.y * rabbits.x, wolves.q ~ γ * wolves.y * rabbits.x]

            ModelingToolkitBase.compose(
                System(eqs, t, [], ps; name = name),
                wolves,
                rabbits,
            )
        end

        function getbyname(sys, name)
            println(name)
            return first(
                [
                    v for v in vcat(unknowns(sys), parameters(sys)) if
                        replace(string(v), "(t)" => "") == name
                ]
            )
        end

        @named ltk_mtk = lotka_volterra_creator()
        simp_ltk_mtk = structural_simplify(ltk_mtk)
        wolves₊δ = getbyname(simp_ltk_mtk, "wolves₊δ")
        rabbits₊α = getbyname(simp_ltk_mtk, "rabbits₊α")
        β = getbyname(simp_ltk_mtk, "β")
        γ = getbyname(simp_ltk_mtk, "γ")
        wolves₊y = getbyname(simp_ltk_mtk, "wolves₊y")
        rabbits₊x = getbyname(simp_ltk_mtk, "rabbits₊x")
        @variables y(t)
        measured_quantities = [y ~ wolves₊y]
        result =
            assess_identifiability(simp_ltk_mtk, measured_quantities = measured_quantities)
        correct = Dict(
            rabbits₊α => :locally,
            γ => :nonidentifiable,
            β => :globally,
            wolves₊δ => :globally,
            rabbits₊x => :nonidentifiable,
            wolves₊y => :globally,
        )
        @test Dict(result) == correct

        #----------------------------------

        @independent_variables t
        @variables x(t), y(t), z(t), w(t)
        @parameters a
        @named sys = System([D(x) ~ a * y], t, [x], [a]; observed = [y ~ z, z ~ x])
        measured_quantities = [w ~ x]
        result = assess_identifiability(sys, measured_quantities = measured_quantities)
        @test result[a] == :globally

        result = find_identifiable_functions(sys, measured_quantities = measured_quantities)
        @test isequal(result, [a])

        #----------------------------------

        # Tensor definition case as reported in
        # https://github.com/SciML/StructuralIdentifiability.jl/issues/178
        @independent_variables t
        @variables x(t)[1:2], y(t)[1:2]
        @parameters k1, k2

        eqs = [D(x[1]) ~ -k1 * x[2], D(x[2]) ~ -k2 * x[1]]

        sys = System(eqs, t, name = :example_vector)
        correct = OrderedDict(x[1] => true, x[2] => true, k1 => true, k2 => true)
        @test assess_local_identifiability(sys, measured_quantities = [x[1], x[2]]) ==
            correct

        # Extension of https://github.com/SciML/StructuralIdentifiability.jl/issues/398
        @parameters k[1:2] a
        @variables (X(t))[1:2] (y(t))[1:2] [output = true]
        eqs = [
            D(X[1]) ~ k[1] - k[2] * X[2],
            D(X[2]) ~ k[1] - k[2] * X[1],
            y[1] ~ X[1] * X[2] + a,
            y[2] ~ X[1] - X[2],
        ]
        @named osys = System(eqs, t)
        correct = OrderedDict(
            X[1] => :locally,
            X[2] => :locally,
            k[1] => :locally,
            k[2] => :globally,
            a => :globally,
        )
        res = assess_identifiability(osys)
        println(res)
        @test res == correct

        #------------------------------------
        # system from the SciML tutorial
        # https://docs.sciml.ai/ModelingToolkit/stable/tutorials/parameter_identifiability/

        @variables x4(t) x5(t) x6(t) x7(t) y1(t) [output = true] y2(t) [output = true]
        @parameters k5, k6, k7, k8, k9, k10
        eqs = [
            D(x4) ~ -k5 * x4 / (k6 + x4),
            D(x5) ~ k5 * x4 / (k6 + x4) - k7 * x5 / (k8 + x5 + x6),
            D(x6) ~ k7 * x5 / (k8 + x5 + x6) - k9 * x6 * (k10 - x6) / k10,
            D(x7) ~ k9 * x6 * (k10 - x6) / k10,
            y1 ~ x4,
            y2 ~ x5,
        ]

        # define the system
        de = System(eqs, t, name = :Biohydrogenation)

        local_id_all = assess_local_identifiability(de, prob_threshold = 0.99)
        @test local_id_all == OrderedDict(
            x4 => true,
            x5 => true,
            x6 => true,
            x7 => false,
            k5 => true,
            k6 => true,
            k7 => true,
            k8 => true,
            k9 => true,
            k10 => true,
        )

        #------------------------------------
        # system from the SciML tutorial
        # https://docs.sciml.ai/ModelingToolkit/stable/tutorials/parameter_identifiability/
        @parameters b, c, α, β, γ, δ, σ
        @variables x1(t) x2(t) x3(t) x4(t) y(t) [output = true] y2(t) [output = true]
        eqs = [
            D(x1) ~ -b * x1 + 1 / (c + x4),
            D(x2) ~ α * x1 - β * x2,
            D(x3) ~ γ * x2 - δ * x3,
            D(x4) ~ σ * x4 * (γ * x2 - δ * x3) / x3,
            y ~ x1 + x2,
            y2 ~ x2,
        ]

        ode = System(eqs, t, name = :GoodwinOsc)

        global_id = assess_identifiability(ode)
        @test global_id == OrderedDict(
            x1 => :globally,
            x2 => :globally,
            x3 => :nonidentifiable,
            x4 => :globally,
            b => :globally,
            c => :globally,
            α => :globally,
            β => :globally,
            γ => :nonidentifiable,
            δ => :globally,
            σ => :globally,
        )
        @test Set(find_identifiable_functions(ode, with_states = true)) ==
            Set([x4, x2, x1, σ, δ, β, α, c, b, x3 / γ])
    end

    @testset "Discrete local identifiability, ModelingToolkit interface" begin
        cases = []

        @independent_variables t
        @parameters α β
        @variables S(t) I(t) R(t) y(t)
        k = ShiftIndex(t)

        eqs = [
            S(k) ~ S(k - 1) - β * S(k - 1) * I(k - 1),
            I(k) ~ I(k - 1) + β * S(k - 1) * I(k - 1) - α * I(k - 1),
            R(k) ~ R(k - 1) + α * I(k - 1),
        ]
        @named sir = System(eqs, t)
        push!(
            cases,
            Dict(
                :dds => sir,
                :res => OrderedDict(S => true, I => true, R => false, α => true, β => true),
                :y => [y ~ I],
                :y2 => [I],
                :known_ic => Array{}[],
                :to_check => Array{}[],
            ),
        )

        @parameters θ
        @variables x(t) y(t)

        eqs = [x(k) ~ θ * x(k - 1)^3]

        @named eqs = System(eqs, t)
        push!(
            cases,
            Dict(
                :dds => eqs,
                :res => OrderedDict(x => true, θ => true),
                :y => [y ~ x],
                :y2 => [x],
                :known_ic => Array{}[],
                :to_check => Array{}[],
            ),
        )

        @parameters θ β
        @variables x1(t) x2(t) y(t)

        eqs = [x1(k) ~ x1(k - 1) + x2(k - 1), x2(k) ~ x2(k - 1) + θ + β]

        @named eqs = System(eqs, t)
        push!(
            cases,
            Dict(
                :dds => eqs,
                :res => OrderedDict(x1 => true, x2 => true, θ => false, β => false),
                :y => [y ~ x1],
                :y2 => [x1],
                :known_ic => Array{}[],
                :to_check => Array{}[],
            ),
        )

        @parameters a b c d
        @variables x1(t) x2(t) y2(t)

        eqs = [
            x1(k) ~ a * x1(k - 1) - b * x1(k - 1) * x2(k - 1),
            x2(k) ~ -c * x2(k - 1) + d * x1(k - 1) * x2(k - 1),
        ]

        @named lv = System(eqs, t)
        push!(
            cases,
            Dict(
                :dds => lv,
                :res => OrderedDict(
                    x1 => true,
                    x2 => false,
                    a => true,
                    b => false,
                    c => true,
                    d => true,
                ),
                :y => [y ~ x1],
                :y2 => [x1],
                :known_ic => Array{}[],
                :to_check => Array{}[],
            ),
        )

        push!(
            cases,
            Dict(
                :dds => lv,
                :res => OrderedDict(b * x2 => true),
                :y => [y ~ x1],
                :y2 => [x1],
                :known_ic => Array{}[],
                :to_check => [b * x2],
            ),
        )

        push!(
            cases,
            Dict(
                :dds => lv,
                :res => OrderedDict(
                    x1 => true,
                    x2 => true,
                    a => true,
                    b => true,
                    c => true,
                    d => true,
                ),
                :y => [y ~ x1, y2 ~ x1 / x2],
                :y2 => [x1, x1 / x2],
                :known_ic => Array{}[],
                :to_check => Array{}[],
            ),
        )

        push!(
            cases,
            Dict(
                :dds => lv,
                :res => OrderedDict(
                    substitute(x1, Dict(t => 0)) => true,
                    substitute(x2, Dict(t => 0)) => true,
                    a => true,
                    b => true,
                    c => true,
                    d => true,
                ),
                :y => [y ~ x1],
                :y2 => [x1],
                :known_ic => [x2],
                :to_check => Array{}[],
            ),
        )

        # Example 1 from https://doi.org/10.1016/j.automatica.2008.03.019
        @parameters theta1 theta2
        @variables x1(t) x2(t) u(t) y(t)

        eqs = [
            x1(k) ~ theta1 * x1(k - 1) + x2(k - 1),
            x2(k) ~ (1 - theta2) * x1(k - 1) + x2(k - 1)^2 + u(k - 1),
        ]

        @named abmd1 = System(eqs, t)
        push!(
            cases,
            Dict(
                :dds => abmd1,
                :res => OrderedDict(x1 => true, x2 => true, theta1 => true, theta2 => true),
                :y => [y ~ x1],
                :y2 => [x1],
                :known_ic => Array{}[],
                :to_check => Array{}[],
            ),
        )

        # Example 2 from https://doi.org/10.1016/j.automatica.2008.03.019
        @parameters theta1 theta2 theta3
        @variables x1(t) x2(t) u(t) y(t) y2(t)

        eqs = [
            x1(k) ~ theta1 * x1(k - 1)^2 + theta2 * x2(k - 1) + u(k - 1),
            x2(k) ~ theta3 * x1(k - 1),
        ]

        @named abmd2 = System(eqs, t)
        push!(
            cases,
            Dict(
                :dds => abmd2,
                :res => OrderedDict(
                    x1 => true,
                    x2 => false,
                    theta1 => true,
                    theta2 => false,
                    theta3 => false,
                ),
                :y => [y ~ x1],
                :y2 => [x1],
                :known_ic => Array{}[],
                :to_check => Array{}[],
            ),
        )
        push!(
            cases,
            Dict(
                :dds => abmd2,
                :res => Dict(
                    x1 => true,
                    x2 => true,
                    theta1 => true,
                    theta2 => true,
                    theta3 => true,
                ),
                :y => [y ~ x1, y2 ~ x2],
                :y2 => [x1, x2],
                :known_ic => Array{}[],
                :to_check => Array{}[],
            ),
        )

        @parameters a b
        @variables x1(t) y(t)

        eqs = [x1(k) ~ x1(k - 1) + a]

        @named kic = System(eqs, t)
        push!(
            cases,
            Dict(
                :dds => kic,
                :res => OrderedDict(x1 => false, a => true, b => false),
                :y => [y ~ x1 + b],
                :y2 => [x1 + b],
                :known_ic => Array{}[],
                :to_check => Array{}[],
            ),
        )
        push!(
            cases,
            Dict(
                :dds => kic,
                :res =>
                    OrderedDict(substitute(x1, Dict(t => 0)) => true, a => true, b => true),
                :y => [y ~ x1 + b],
                :y2 => [x1 + b],
                :known_ic => [x1],
                :to_check => Array{}[],
            ),
        )

        for c in cases
            @test assess_local_identifiability(
                c[:dds];
                measured_quantities = c[:y],
                known_ic = c[:known_ic],
                funcs_to_check = c[:to_check],
            ) == c[:res]
            @test assess_local_identifiability(
                c[:dds];
                measured_quantities = c[:y2],
                known_ic = c[:known_ic],
                funcs_to_check = c[:to_check],
            ) == c[:res]
        end
    end

    @testset "Exporting ModelingToolkit Model to SI Model" begin

        # Creates MTK model and assesses its identifiability.
        @independent_variables t
        @parameters r1, r2, c1, c2, beta1, beta2, chi1, chi2
        @variables x1(t), x2(t), y(t), u(t)
        D = Differential(t)
        eqs = [
            D(x1) ~ r1 * x1 * (1 - c1 * x1) + beta1 * x1 * x2 / (chi1 + x2) + u,
            D(x2) ~ r2 * x2 * (1 - c2 * x2) + beta2 * x1 * x2 / (chi2 + x1),
        ]
        measured_quantities = [y ~ x1]
        ode_mtk = System(eqs, t, name = :mutualist)

        global_id_1 =
            assess_identifiability(ode_mtk, measured_quantities = measured_quantities)
        local_id_1 =
            assess_local_identifiability(ode_mtk, measured_quantities = measured_quantities)
        ifs_1 =
            find_identifiable_functions(ode_mtk, measured_quantities = measured_quantities)

        # Converts mtk model to si model, and assesses its identifiability.
        si_model, _ = mtk_to_si(ode_mtk, measured_quantities)
        global_id_2 = assess_identifiability(si_model)
        local_id_2 = assess_local_identifiability(si_model)
        ifs_2 = find_identifiable_functions(si_model)

        # Converts the output dicts from StructuralIdentifiability functions from "weird symbol => stuff" to "symbol => stuff" (the output have some strange meta data which prevents equality checks, this enables this).
        # Structural identifiability also provides variables like x (rather than x(t)). This is a bug, but we have to convert to make it work (now just remove any (t) to make them all equal).
        function sym_dict(dict_in)
            dict_out = Dict{Symbol, Any}()
            for key in keys(dict_in)
                sym_key = Symbol(key)
                sym_key = Symbol(replace(String(sym_key), "(t)" => ""))
                dict_out[sym_key] = dict_in[key]
            end
            return dict_out
        end

        println(sym_dict(local_id_1))
        println(sym_dict(local_id_2))
        # Checks that the two approaches yields the same result
        @test sym_dict(local_id_1) == sym_dict(local_id_2)
        @test sym_dict(global_id_1) == sym_dict(global_id_2)
        @test length(ifs_1) == length(ifs_2)

        # Does not take `observables` as outputs if there is something else
        @variables X(t) μ₁(t) yX(t) [output = true]
        @parameters k₁ k₂ μ₁max μ₂max
        eqs = [μ₁ ~ k₁ + μ₁max, D(X) ~ (μ₁ + μ₂max) * X, yX ~ X]

        ode = System(eqs, t, name = :output_definition_case)

        id_res = assess_identifiability(ode)
        @test 1 == count(v -> v == :globally, values(id_res))
    end

    @testset "Identifiability of MTK models with known generic initial conditions" begin
        cases = []

        @independent_variables t
        @parameters a, b, c, d
        @variables x1(t), x2(t)
        D = Differential(t)
        x1_0 = substitute(x1, Dict(t => 0))
        x2_0 = substitute(x2, Dict(t => 0))
        eqs = [D(x1) ~ a * x1 - b * x1 * x2, D(x2) ~ -c * x2 + d * x1 * x2]
        ode_mtk = System(eqs, t, name = :lv)
        push!(
            cases,
            Dict(
                :ode => ode_mtk,
                :measured => [x1],
                :known => [x2],
                :to_check => [],
                :correct_funcs => [a, b, c, d, x1_0, x2_0],
                :correct_ident =>
                    OrderedDict(x => :globally for x in [x1_0, x2_0, a, b, c, d]),
            ),
        )

        @parameters c
        @variables x3(t)
        x3_0 = substitute(x3, Dict(t => 0))
        eqs = [D(x1) ~ a + x2 + x3, D(x2) ~ b^2 + c, D(x3) ~ -c]
        ode_mtk = System(eqs, t, name = :ex2)

        push!(
            cases,
            Dict(
                :ode => ode_mtk,
                :measured => [x1],
                :known => [x2, x3],
                :to_check => [],
                :correct_funcs => [a, b^2, x1_0, x2_0, x3_0],
                :correct_ident => OrderedDict(
                    x1_0 => :globally,
                    x2_0 => :globally,
                    x3_0 => :globally,
                    a => :globally,
                    b => :locally,
                    c => :nonidentifiable,
                ),
            ),
        )

        push!(
            cases,
            Dict(
                :ode => ode_mtk,
                :measured => [x1],
                :known => [x2, x3],
                :to_check => [b^2, x2 * c],
                :correct_funcs => [a, b^2, x1_0, x2_0, x3_0],
                :correct_ident =>
                    OrderedDict(b^2 => :globally, x2_0 * c => :nonidentifiable),
            ),
        )

        eqs = [D(x1) ~ a * x1, D(x2) ~ x1 + 1 / x1]
        ode_mtk = System(eqs, t, name = :ex3)
        push!(
            cases,
            Dict(
                :ode => ode_mtk,
                :measured => [x2],
                :known => [x1],
                :to_check => [],
                :correct_funcs => [a, x1_0, x2_0],
                :correct_ident =>
                    OrderedDict(x1_0 => :globally, x2_0 => :globally, a => :globally),
            ),
        )

        @parameters alpha, beta, gama, delta, sigma
        @variables x4(t)
        x4_0 = substitute(x4, Dict(t => 0))
        eqs = [
            D(x1) ~ -b * x1 + 1 / (c + x4),
            D(x2) ~ alpha * x1 - beta * x2,
            D(x3) ~ gama * x2 - delta * x3,
            D(x4) ~ sigma * x4 * (gama * x2 - delta * x3) / x3,
        ]
        ode_mtk = System(eqs, t, name = :goodwin)
        push!(
            cases,
            Dict(
                :ode => ode_mtk,
                :measured => [x1],
                :known => [x2, x3],
                :to_check => [alpha, alpha * gama],
                :correct_funcs => [
                    sigma,
                    c,
                    b,
                    x4_0,
                    x3_0,
                    x2_0,
                    x1_0,
                    beta * delta,
                    alpha * gama,
                    beta + delta,
                    -delta * x3_0 + gama * x2_0,
                ],
                :correct_ident => OrderedDict(alpha => :locally, alpha * gama => :globally),
            ),
        )

        for case in cases
            ode = case[:ode]
            y = case[:measured]
            known = case[:known]
            result_funcs =
                find_identifiable_functions(ode, known_ic = known, measured_quantities = y)
            correct_funcs = @test Set(result_funcs) == Set(case[:correct_funcs])

            result_ident = assess_identifiability(
                ode,
                known_ic = known,
                measured_quantities = y,
                funcs_to_check = case[:to_check],
            )
            @test case[:correct_ident] == result_ident
        end
    end
end
