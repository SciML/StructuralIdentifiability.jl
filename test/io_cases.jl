@testset "Checking io-equations: single output" begin
    test_cases = []

    # 2-compartiment model
    push!(test_cases,
          Dict("ode" => @ODEmodel(x0'(t)=-(a01 + a21) * x0(t) + a12 * x1(t) + u(t),
                                  x1'(t)=a21 * x0(t) - a12 * x1(t),
                                  y(t)=x0(t)),
               "correct" => "y_2 + y_1*a01 + y_1*a21 + y_1*a12 + y_0*a01*a12 - u_0*a12 - u_1"))

    #---------------------------------------
    # Chen-Lee model
    correct = "9*y1_3^2*y1_0^2 - 18*y1_3*y1_2*y1_1*y1_0 - 18*y1_3*y1_2*y1_0^2*a - 36*y1_3*y1_2*y1_0^2*b - 36*y1_3*y1_2*y1_0^2*c + 18*y1_3*y1_1^2*y1_0*a + 18*y1_3*y1_1^2*y1_0*b + 18*y1_3*y1_1^2*y1_0*c - 24*y1_3*y1_1*y1_0^4 + 18*y1_3*y1_1*y1_0^2*a*b + 18*y1_3*y1_1*y1_0^2*a*c + 18*y1_3*y1_1*y1_0^2*b^2 + 36*y1_3*y1_1*y1_0^2*b*c + 18*y1_3*y1_1*y1_0^2*c^2 + 24*y1_3*y1_0^5*a - 18*y1_3*y1_0^3*a*b^2 - 36*y1_3*y1_0^3*a*b*c - 18*y1_3*y1_0^3*a*c^2 + 9*y1_2^2*y1_1^2 + 18*y1_2^2*y1_1*y1_0*a + 36*y1_2^2*y1_1*y1_0*b + 36*y1_2^2*y1_1*y1_0*c + 9*y1_2^2*y1_0^2*a^2 + 36*y1_2^2*y1_0^2*a*b + 36*y1_2^2*y1_0^2*a*c + 27*y1_2^2*y1_0^2*b^2 + 90*y1_2^2*y1_0^2*b*c + 27*y1_2^2*y1_0^2*c^2 - 18*y1_2*y1_1^3*a - 18*y1_2*y1_1^3*b - 18*y1_2*y1_1^3*c + 24*y1_2*y1_1^2*y1_0^3 - 18*y1_2*y1_1^2*y1_0*a^2 - 72*y1_2*y1_1^2*y1_0*a*b - 72*y1_2*y1_1^2*y1_0*a*c - 54*y1_2*y1_1^2*y1_0*b^2 - 108*y1_2*y1_1^2*y1_0*b*c - 54*y1_2*y1_1^2*y1_0*c^2 + 48*y1_2*y1_1*y1_0^4*b + 48*y1_2*y1_1*y1_0^4*c - 18*y1_2*y1_1*y1_0^2*a^2*b - 18*y1_2*y1_1*y1_0^2*a^2*c - 18*y1_2*y1_1*y1_0^2*a*b^2 - 108*y1_2*y1_1*y1_0^2*a*b*c - 18*y1_2*y1_1*y1_0^2*a*c^2 - 18*y1_2*y1_1*y1_0^2*b^3 - 126*y1_2*y1_1*y1_0^2*b^2*c - 126*y1_2*y1_1*y1_0^2*b*c^2 - 18*y1_2*y1_1*y1_0^2*c^3 - 24*y1_2*y1_0^5*a^2 - 48*y1_2*y1_0^5*a*b - 48*y1_2*y1_0^5*a*c + 18*y1_2*y1_0^3*a^2*b^2 + 36*y1_2*y1_0^3*a^2*b*c + 18*y1_2*y1_0^3*a^2*c^2 + 18*y1_2*y1_0^3*a*b^3 + 126*y1_2*y1_0^3*a*b^2*c + 126*y1_2*y1_0^3*a*b*c^2 + 18*y1_2*y1_0^3*a*c^3 + 9*y1_1^4*a^2 + 18*y1_1^4*a*b + 18*y1_1^4*a*c + 9*y1_1^4*b^2 + 18*y1_1^4*b*c + 9*y1_1^4*c^2 - 24*y1_1^3*y1_0^3*a - 24*y1_1^3*y1_0^3*b - 24*y1_1^3*y1_0^3*c + 18*y1_1^3*y1_0*a^2*b + 18*y1_1^3*y1_0*a^2*c + 36*y1_1^3*y1_0*a*b^2 + 72*y1_1^3*y1_0*a*b*c + 36*y1_1^3*y1_0*a*c^2 + 18*y1_1^3*y1_0*b^3 + 54*y1_1^3*y1_0*b^2*c + 54*y1_1^3*y1_0*b*c^2 + 18*y1_1^3*y1_0*c^3 + 16*y1_1^2*y1_0^6 + 24*y1_1^2*y1_0^4*a^2 - 12*y1_1^2*y1_0^4*b^2 - 72*y1_1^2*y1_0^4*b*c - 12*y1_1^2*y1_0^4*c^2 - 18*y1_1^2*y1_0^2*a^2*b^2 - 18*y1_1^2*y1_0^2*a^2*c^2 - 18*y1_1^2*y1_0^2*a*b^3 + 18*y1_1^2*y1_0^2*a*b^2*c + 18*y1_1^2*y1_0^2*a*b*c^2 - 18*y1_1^2*y1_0^2*a*c^3 + 36*y1_1^2*y1_0^2*b^3*c + 72*y1_1^2*y1_0^2*b^2*c^2 + 36*y1_1^2*y1_0^2*b*c^3 - 32*y1_1*y1_0^7*a + 24*y1_1*y1_0^5*a^2*b + 24*y1_1*y1_0^5*a^2*c + 24*y1_1*y1_0^5*a*b^2 + 144*y1_1*y1_0^5*a*b*c + 24*y1_1*y1_0^5*a*c^2 - 72*y1_1*y1_0^3*a^2*b^2*c - 72*y1_1*y1_0^3*a^2*b*c^2 - 72*y1_1*y1_0^3*a*b^3*c - 144*y1_1*y1_0^3*a*b^2*c^2 - 72*y1_1*y1_0^3*a*b*c^3 + 16*y1_0^8*a^2 - 12*y1_0^6*a^2*b^2 - 72*y1_0^6*a^2*b*c - 12*y1_0^6*a^2*c^2 + 36*y1_0^4*a^2*b^3*c + 72*y1_0^4*a^2*b^2*c^2 + 36*y1_0^4*a^2*b*c^3"

    push!(test_cases,
          Dict("ode" => @ODEmodel(x0'(t)=a * x0(t) - x1(t) * x2(t),
                                  x1'(t)=b * x1(t) + x0(t) * x2(t),
                                  x2'(t)=c * x2(t) + 1 // 3 * x0(t) * x1(t),
                                  y1(t)=x0(t)),
               "correct" => correct))

    #---------------------------------------

    # predator-prey model
    correct = "y1_2*y1_0 - y1_1^2 - y1_1*y1_0^2*d + y1_1*y1_0*c + y1_0^3*a*d - y1_0^2*a*c"

    push!(test_cases,
          Dict("ode" => @ODEmodel(x0'(t)=a * x0(t) - b * x0(t) * x1(t),
                                  x1'(t)=-c * x1(t) + d * x0(t) * x1(t),
                                  y1(t)=x0(t)),
               "correct" => correct))

    #---------------------------------------

    for case in test_cases
        ode = case["ode"]
        io_eq = collect(values(find_ioequations(ode)))[1]
        gen_symbols = symbols(parent(io_eq))
        correct_expr = Meta.parse(case["correct"])
        for (i, s) in enumerate(gen_symbols)
            correct_expr = StructuralIdentifiability.MacroTools.postwalk(ex -> ex == s ?
                                                                               :(polyvars[$i]) :
                                                                               ex,
                                                                         correct_expr)
        end
        # dirty hack to eval within the local scope
        correct_function = @eval function (eq)
            polyvars = gens(parent(eq))
            $correct_expr
        end
        correct = correct_function(io_eq)
        divisibility, remainder = divides(io_eq, correct)
        @test divisibility
        @test total_degree(remainder) == 0
    end
end
