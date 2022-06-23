@testset "Assessing local identifiability" begin
    test_cases = []

    # 2-compartiment model

    ode = @ODEmodel(x0'(t) = -(a01 + a21) * x0(t) + a12 * x1(t),
                    x1'(t) = a21 * x0(t) - a12 * x1(t),
                    y(t) = x0(t))
    funcs_to_test = [a01, a21, a12, a01 * a12, a01 + a12 + a21, (a01 + a12 + a21) // (a01 * a12),
                     x0, x1]
    correct = [false, false, false, true, true, true, true, false]
    push!(test_cases,
          Dict(:ode => ode,
               :funcs => funcs_to_test,
               :correct => Dict(funcs_to_test .=> correct)))

    #--------------------------------------------------------------------------

    ode = @ODEmodel(x0'(t) = a * x0(t) - b * x0(t) * x1(t) + u(t),
                    x1'(t) = c * x1(t) + d * x0(t) * x1(t),
                    y(t) = x0(t))
    funcs_to_test = [a, b, c, d, x0, x1]
    correct = [true, false, true, true, true, false]
    push!(test_cases,
          Dict(:ode => ode,
               :funcs => funcs_to_test,
               :correct => Dict(funcs_to_test .=> correct)))

    #--------------------------------------------------------------------------

    ode = @ODEmodel(S'(t) = mu - bi * S(t) * I(t) - bw * S(t) * W(t) - mu * S(t) + a * R(t),
                    I'(t) = bw * S(t) * W(t) + bi * S(t) * I(t) - (gam + mu) * I(t),
                    W'(t) = xi * (I(t) - W(t)),
                    R'(t) = gam * I(t) - (mu + a) * R(t),
                    y(t) = k * I(t))
    funcs_to_test = [mu, bi, bw, a, xi, gam, mu, gam + mu, k, S, I, W, R]
    correct = [true for _ in funcs_to_test]
    push!(test_cases,
          Dict(:ode => ode,
               :funcs => funcs_to_test,
               :correct => Dict(funcs_to_test .=> correct)))

    #--------------------------------------------------------------------------

    ode = @ODEmodel(x1'(t) = -b * x1(t) + 1 / (c + x4(t)),
                    x2'(t) = alpha * x1(t) - beta * x2(t),
                    x3'(t) = gama * x2(t) - delta * x3(t),
                    x4'(t) = sigma * x4(t) * (gama * x2(t) - delta * x3(t)) / x3(t),
                    y(t) = x1(t))
    funcs_to_test = [b, c, alpha, beta, delta, gama, beta + delta, beta * delta]
    correct = Dict([b => true, c => true, alpha => false, beta => true, delta => true,
                    gama => false, beta + delta => true, beta * delta => true])
    push!(test_cases, Dict(:ode => ode,
                           :funcs => funcs_to_test,
                           :correct => correct))

    #--------------------------------------------------------------------------

    for case in test_cases
        result = assess_local_identifiability(case[:ode], case[:funcs])
        @test result == case[:correct]
    end
end
