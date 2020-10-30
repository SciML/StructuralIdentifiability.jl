using MacroTools

@testset "Checking io-equations: single output" begin
    test_cases = []

    P = fmpq_mpoly

    # 2-compartiment model
    var_names = [
        "x_0", "x_1", "u",
        "a_01", "a_21", "a_12"
    ]
    R, (x_0, x_1, u, a_01, a_21, a_12) = PolynomialRing(QQ, var_names)
    f = [-(a_01 + a_21) * x_0 + a_12 * x_1 + u, a_21 * x_0 - a_12 * x_1]
    g = x_0

    ode = ODE{P}(Dict{P, Union{P, Generic.Frac{P}}}(x_0 => f[1], x_1 => f[2]), [u])
    push!(test_cases, Dict(
        "ode" => ode,
        "output" => g,
        "correct" => "y1_2 + y1_1*a_01 + y1_1*a_21 + y1_1*a_12 + y1_0*a_01*a_12 - u*a_12 - u_1"
    ))

    #---------------------------------------
    # Chen-Lee model
    var_names = [
        "x_0", "x_1", "x_2",
        "a", "b", "c"
    ]
    R, (x_0, x_1, x_2, a, b, c) = PolynomialRing(QQ, var_names)
    f = [
        a * x_0 - x_1 * x_2, 
        b*x_1 + x_0 * x_2, 
        c * x_2 + 1//3 * x_0 * x_1
    ]
    g = x_0

    ode = ODE{P}(Dict{P, Union{P, Generic.Frac{P}}}(x_0 => f[1], x_1 => f[2], x_2 => f[3]), Array{P, 1}())

    correct = "9*y1_3^2*y1_0^2 - 18*y1_3*y1_2*y1_1*y1_0 - 18*y1_3*y1_2*y1_0^2*a - 36*y1_3*y1_2*y1_0^2*b - 36*y1_3*y1_2*y1_0^2*c + 18*y1_3*y1_1^2*y1_0*a + 18*y1_3*y1_1^2*y1_0*b + 18*y1_3*y1_1^2*y1_0*c - 24*y1_3*y1_1*y1_0^4 + 18*y1_3*y1_1*y1_0^2*a*b + 18*y1_3*y1_1*y1_0^2*a*c + 18*y1_3*y1_1*y1_0^2*b^2 + 36*y1_3*y1_1*y1_0^2*b*c + 18*y1_3*y1_1*y1_0^2*c^2 + 24*y1_3*y1_0^5*a - 18*y1_3*y1_0^3*a*b^2 - 36*y1_3*y1_0^3*a*b*c - 18*y1_3*y1_0^3*a*c^2 + 9*y1_2^2*y1_1^2 + 18*y1_2^2*y1_1*y1_0*a + 36*y1_2^2*y1_1*y1_0*b + 36*y1_2^2*y1_1*y1_0*c + 9*y1_2^2*y1_0^2*a^2 + 36*y1_2^2*y1_0^2*a*b + 36*y1_2^2*y1_0^2*a*c + 27*y1_2^2*y1_0^2*b^2 + 90*y1_2^2*y1_0^2*b*c + 27*y1_2^2*y1_0^2*c^2 - 18*y1_2*y1_1^3*a - 18*y1_2*y1_1^3*b - 18*y1_2*y1_1^3*c + 24*y1_2*y1_1^2*y1_0^3 - 18*y1_2*y1_1^2*y1_0*a^2 - 72*y1_2*y1_1^2*y1_0*a*b - 72*y1_2*y1_1^2*y1_0*a*c - 54*y1_2*y1_1^2*y1_0*b^2 - 108*y1_2*y1_1^2*y1_0*b*c - 54*y1_2*y1_1^2*y1_0*c^2 + 48*y1_2*y1_1*y1_0^4*b + 48*y1_2*y1_1*y1_0^4*c - 18*y1_2*y1_1*y1_0^2*a^2*b - 18*y1_2*y1_1*y1_0^2*a^2*c - 18*y1_2*y1_1*y1_0^2*a*b^2 - 108*y1_2*y1_1*y1_0^2*a*b*c - 18*y1_2*y1_1*y1_0^2*a*c^2 - 18*y1_2*y1_1*y1_0^2*b^3 - 126*y1_2*y1_1*y1_0^2*b^2*c - 126*y1_2*y1_1*y1_0^2*b*c^2 - 18*y1_2*y1_1*y1_0^2*c^3 - 24*y1_2*y1_0^5*a^2 - 48*y1_2*y1_0^5*a*b - 48*y1_2*y1_0^5*a*c + 18*y1_2*y1_0^3*a^2*b^2 + 36*y1_2*y1_0^3*a^2*b*c + 18*y1_2*y1_0^3*a^2*c^2 + 18*y1_2*y1_0^3*a*b^3 + 126*y1_2*y1_0^3*a*b^2*c + 126*y1_2*y1_0^3*a*b*c^2 + 18*y1_2*y1_0^3*a*c^3 + 9*y1_1^4*a^2 + 18*y1_1^4*a*b + 18*y1_1^4*a*c + 9*y1_1^4*b^2 + 18*y1_1^4*b*c + 9*y1_1^4*c^2 - 24*y1_1^3*y1_0^3*a - 24*y1_1^3*y1_0^3*b - 24*y1_1^3*y1_0^3*c + 18*y1_1^3*y1_0*a^2*b + 18*y1_1^3*y1_0*a^2*c + 36*y1_1^3*y1_0*a*b^2 + 72*y1_1^3*y1_0*a*b*c + 36*y1_1^3*y1_0*a*c^2 + 18*y1_1^3*y1_0*b^3 + 54*y1_1^3*y1_0*b^2*c + 54*y1_1^3*y1_0*b*c^2 + 18*y1_1^3*y1_0*c^3 + 16*y1_1^2*y1_0^6 + 24*y1_1^2*y1_0^4*a^2 - 12*y1_1^2*y1_0^4*b^2 - 72*y1_1^2*y1_0^4*b*c - 12*y1_1^2*y1_0^4*c^2 - 18*y1_1^2*y1_0^2*a^2*b^2 - 18*y1_1^2*y1_0^2*a^2*c^2 - 18*y1_1^2*y1_0^2*a*b^3 + 18*y1_1^2*y1_0^2*a*b^2*c + 18*y1_1^2*y1_0^2*a*b*c^2 - 18*y1_1^2*y1_0^2*a*c^3 + 36*y1_1^2*y1_0^2*b^3*c + 72*y1_1^2*y1_0^2*b^2*c^2 + 36*y1_1^2*y1_0^2*b*c^3 - 32*y1_1*y1_0^7*a + 24*y1_1*y1_0^5*a^2*b + 24*y1_1*y1_0^5*a^2*c + 24*y1_1*y1_0^5*a*b^2 + 144*y1_1*y1_0^5*a*b*c + 24*y1_1*y1_0^5*a*c^2 - 72*y1_1*y1_0^3*a^2*b^2*c - 72*y1_1*y1_0^3*a^2*b*c^2 - 72*y1_1*y1_0^3*a*b^3*c - 144*y1_1*y1_0^3*a*b^2*c^2 - 72*y1_1*y1_0^3*a*b*c^3 + 16*y1_0^8*a^2 - 12*y1_0^6*a^2*b^2 - 72*y1_0^6*a^2*b*c - 12*y1_0^6*a^2*c^2 + 36*y1_0^4*a^2*b^3*c + 72*y1_0^4*a^2*b^2*c^2 + 36*y1_0^4*a^2*b*c^3"

    push!(test_cases, Dict(
        "ode" => ode,
        "output" => g,
        "correct" => correct
    ))

    #---------------------------------------

    # predator-prey model
    var_names = [
        "x_0", "x_1",
        "a", "b", "c", "d"
    ]
    R, (x_0, x_1, a, b, c, d) = PolynomialRing(QQ, var_names)
    f = [a * x_0 - b * x_0 * x_1, - c * x_1 + d * x_0 * x_1]
    g = x_0
    ode = ODE{P}(Dict{P, Union{P, Generic.Frac{P}}}(x_0 => f[1], x_1 => f[2]), Array{P, 1}())

    correct = "y1_2*y1_0 - y1_1^2 - y1_1*y1_0^2*d + y1_1*y1_0*c + y1_0^3*a*d - y1_0^2*a*c"
 
    push!(test_cases, Dict(
        "ode" => ode,
        "output" => g,
        "correct" => correct
    ))

    #---------------------------------------

    for case in test_cases
        ode = case["ode"]
        g = case["output"]
        io_eq = collect(values(find_ioequations(ode, [g])))[1]
        gen_symbols = symbols(parent(io_eq))
        correct_expr = Meta.parse(case["correct"])
        for (i, s) in enumerate(gen_symbols)
            correct_expr = MacroTools.postwalk(ex -> ex == s ? :(polyvars[$i]) : ex, correct_expr)
        end
        # dirty hack to eval within the local scope
        correct_function = @eval function(eq)
            polyvars = gens(parent(eq))
            $correct_expr
        end
        correct = correct_function(io_eq)
        divisibility, remainder = divides(io_eq, correct)
        @test divisibility
        @test total_degree(remainder) == 0
    end
end
