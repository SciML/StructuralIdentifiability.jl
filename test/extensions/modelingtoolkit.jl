if GROUP == "All" || GROUP == "ModelingToolkitExt"
    @testset "Check identifiability of `ODESystem` object" begin
        using ModelingToolkit
        using ModelingToolkit: parameters
        using Symbolics

        @parameters a01 a21 a12
        @variables t x0(t) x1(t) y1(t) [output = true]
        D = Differential(t)

        eqs =
            [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1, y1 ~ x0]
        de = ODESystem(eqs, t, name = :Test)

        correct = OrderedDict(
            a01 => :nonidentifiable,
            a21 => :nonidentifiable,
            a12 => :nonidentifiable,
            x0 => :globally,
            x1 => :nonidentifiable,
        )

        @test isequal(
            correct,
            assess_identifiability(de; measured_quantities = [y1 ~ x0]),
        )
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

        correct = [k01, c * k_m, V_m * c]
        result = find_identifiable_functions(de)
        @test isequal(Set(correct), Set(result))

        correct = [k01, c * x, k_m * c, V_m * c]
        result = find_identifiable_functions(de, with_states = true)
        @test isequal(Set(correct), Set(result))

        # --------------------------------------------------------------------------
        @parameters a01 a21 a12
        @variables t x0(t) x1(t) y1(t) [output = true]
        D = Differential(t)

        eqs =
            [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1, y1 ~ x0]
        de = ODESystem(eqs, t, name = :Test)

        correct = OrderedDict(
            a01 => :nonidentifiable,
            a21 => :nonidentifiable,
            a12 => :nonidentifiable,
            x0 => :globally,
            x1 => :nonidentifiable,
        )

        @test isequal(correct, assess_identifiability(de))

        # --------------------------------------------------------------------------

        @parameters a01 a21 a12
        @variables t x0(t) x1(t) y1(t) [output = true]
        D = Differential(t)

        eqs =
            [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1, y1 ~ x0]
        de = ODESystem(eqs, t, name = :Test)
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
            assess_identifiability(de; funcs_to_check = funcs_to_check),
        )

        # --------------------------------------------------------------------------

        @parameters a01 a21 a12
        @variables t x0(t) x1(t) y1(t)
        D = Differential(t)

        eqs = [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1]
        measured_quantities = [y1 ~ x0]
        de = ODESystem(eqs, t, name = :Test)
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
        result =
            find_identifiable_functions(de, measured_quantities = measured_quantities)
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
            all(
                values(
                    assess_local_identifiability(de; measured_quantities = [y ~ k * I]),
                ),
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
        result =
            find_identifiable_functions(de, measured_quantities = measured_quantities)
        correct = [b, a, c^2]
        @test isequal(Set(result), Set(correct))

        # ----------
        @parameters a b
        @variables t c(t) x1(t) x2(t) y1(t) y2(t)
        D = Differential(t)

        eqs = [
            D(x1) ~ -a * x1 + x2 * b / (x1 + b / (c^2 - x2)),
            D(x2) ~ x2 * c^2 + x1,
            D(c) ~ 0,
        ]
        de = ODESystem(eqs, t, name = :Test)
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

            eqs = [
                rabbits.z ~ -β * wolves.y * rabbits.x,
                wolves.q ~ γ * wolves.y * rabbits.x,
            ]

            ModelingToolkit.compose(
                ODESystem(eqs, t, [], ps; name = name),
                wolves,
                rabbits,
            )
        end

        function getbyname(sys, name)
            println(name)
            return first([
                v for v in vcat(states(sys), parameters(sys)) if
                replace(string(v), "(t)" => "") == name
            ])
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
        result = assess_identifiability(
            simp_ltk_mtk,
            measured_quantities = measured_quantities,
        )
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

        @variables t, x(t), y(t), z(t), w(t)
        @parameters a
        @named sys = ODESystem([D(x) ~ a * y], t, [x], [a]; observed = [y ~ z, z ~ x])
        measured_quantities = [w ~ x]
        result = assess_identifiability(sys, measured_quantities = measured_quantities)
        @test result[a] == :globally

        result =
            find_identifiable_functions(sys, measured_quantities = measured_quantities)
        @test isequal(result, [a])

        #----------------------------------

        # Tensor definition case as reported in
        # https://github.com/SciML/StructuralIdentifiability.jl/issues/178
        @variables t, x(t)[1:2], y(t)[1:2]
        @parameters k1, k2

        eqs = [D(x[1]) ~ -k1 * x[2], D(x[2]) ~ -k2 * x[1]]

        sys = ODESystem(eqs, t, name = :example_vector)
        correct = OrderedDict(x[1] => true, x[2] => true, k1 => true, k2 => true)
        @test assess_local_identifiability(sys, measured_quantities = [x[1], x[2]]) ==
              correct
    end

    @testset "Discrete local identifiability, ModelingToolkit interface" begin
        cases = []

        @parameters α β
        @variables t S(t) I(t) R(t) y(t)
        D = Difference(t; dt = 1.0)

        eqs = [D(S) ~ -β * S * I, D(I) ~ β * S * I - α * I, D(R) ~ α * I]
        @named sir = DiscreteSystem(eqs)
        push!(
            cases,
            Dict(
                :dds => sir,
                :res =>
                    OrderedDict(S => true, I => true, R => false, α => true, β => true),
                :y => [y ~ I],
                :y2 => [I],
                :known_ic => Array{}[],
                :to_check => Array{}[],
            ),
        )

        @parameters θ
        @variables t x(t) y(t)
        D = Difference(t; dt = 1.0)

        eqs = [D(x) ~ θ * x^3 - x]

        @named eqs = DiscreteSystem(eqs)
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
        @variables t x1(t) x2(t) y(t)
        D = Difference(t; dt = 1.0)

        eqs = [D(x1) ~ x1 + x2, D(x2) ~ θ + β]

        @named eqs = DiscreteSystem(eqs)
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
        @variables t x1(t) x2(t) u(t) y2(t)
        D = Difference(t; dt = 1.0)

        eqs = [D(x1) ~ a * x1 - b * x1 * x2 + u, D(x2) ~ -c * x2 + d * x1 * x2]

        @named lv = DiscreteSystem(eqs)
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
        @variables t x1(t) x2(t) u(t) y(t)
        D = Difference(t; dt = 1.0)

        eqs = [D(x1) ~ theta1 * x1 + x2, D(x2) ~ (1 - theta2) * x1 + x2^2 + u - x2]

        @named abmd1 = DiscreteSystem(eqs)
        push!(
            cases,
            Dict(
                :dds => abmd1,
                :res =>
                    OrderedDict(x1 => true, x2 => true, theta1 => true, theta2 => true),
                :y => [y ~ x1],
                :y2 => [x1],
                :known_ic => Array{}[],
                :to_check => Array{}[],
            ),
        )

        # Example 2 from https://doi.org/10.1016/j.automatica.2008.03.019
        @parameters theta1 theta2 theta3
        @variables t x1(t) x2(t) u(t) y(t) y2(t)
        D = Difference(t; dt = 1.0)

        eqs = [D(x1) ~ theta1 * x1^2 + theta2 * x2 + u - x1, D(x2) ~ theta3 * x1 - x2]

        @named abmd2 = DiscreteSystem(eqs)
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
        @variables t x1(t) y(t)
        D = Difference(t; dt = 1.0)

        eqs = [D(x1) ~ a]

        @named kic = DiscreteSystem(eqs)
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
                :res => OrderedDict(
                    substitute(x1, Dict(t => 0)) => true,
                    a => true,
                    b => true,
                ),
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
end
