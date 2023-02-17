@testset "Assessing identifiability" begin
    test_cases = []

    # 2-compartiment model

    ode = @ODEmodel(x0'(t)=-(a01 + a21) * x0(t) + a12 * x1(t),
                    x1'(t)=a21 * x0(t) - a12 * x1(t),
                    y(t)=x0(t))
    funcs_to_test = [
        a01,
        a21,
        a12,
        a01 * a12,
        a01 + a12 + a21,
        (a01 + a12 + a21) // (a01 * a12),
    ]
    correct = [
        :nonidentifiable,
        :nonidentifiable,
        :nonidentifiable,
        :globally,
        :globally,
        :globally,
    ]
    push!(test_cases,
          Dict(:ode => ode,
               :funcs => funcs_to_test,
               :correct => Dict(funcs_to_test .=> correct)))

    #--------------------------------------------------------------------------

    ode = @ODEmodel(x0'(t)=a * x0(t) - b * x0(t) * x1(t) + u(t),
                    x1'(t)=c * x1(t) + d * x0(t) * x1(t),
                    y(t)=x0(t))
    funcs_to_test = [a, b, c, d]
    correct = [:globally, :nonidentifiable, :globally, :globally]
    push!(test_cases,
          Dict(:ode => ode,
               :funcs => funcs_to_test,
               :correct => Dict(funcs_to_test .=> correct)))

    #--------------------------------------------------------------------------

    ode = @ODEmodel(S'(t)=mu - bi * S(t) * I(t) - bw * S(t) * W(t) - mu * S(t) + a * R(t),
                    I'(t)=bw * S(t) * W(t) + bi * S(t) * I(t) - (gam + mu) * I(t),
                    W'(t)=xi * (I(t) - W(t)),
                    R'(t)=gam * I(t) - (mu + a) * R(t),
                    y(t)=k * I(t))
    funcs_to_test = [mu, bi, bw, a, xi, gam, mu, gam + mu, k]
    correct = [:globally for _ in funcs_to_test]
    push!(test_cases,
          Dict(:ode => ode,
               :funcs => funcs_to_test,
               :correct => Dict(funcs_to_test .=> correct)))

    #--------------------------------------------------------------------------

    ode = @ODEmodel(x1'(t)=-b * x1(t) + 1 / (c + x4(t)),
                    x2'(t)=alpha * x1(t) - beta * x2(t),
                    x3'(t)=gama * x2(t) - delta * x3(t),
                    x4'(t)=sigma * x4(t) * (gama * x2(t) - delta * x3(t)) / x3(t),
                    y(t)=x1(t))
    funcs_to_test = [b, c, alpha, beta, delta, gama, beta + delta, beta * delta]
    correct = [
        :globally,
        :globally,
        :nonidentifiable,
        :locally,
        :locally,
        :nonidentifiable,
        :globally,
        :globally,
    ]
    push!(test_cases,
          Dict(:ode => ode,
               :funcs => funcs_to_test,
               :correct => Dict(funcs_to_test .=> correct)))

    #--------------------------------------------------------------------------

    ode = @ODEmodel(x1'(t)=(1 + x1(t)^2) // 2,
                    x2'(t)=(1 - x1(t)^2) // (1 + x1(t)^2),
                    y1(t)=2 * x1(t) // (b * (1 + x1(t)^2)),
                    y2(t)=x2(t))
    funcs_to_test = [b]
    correct = [:globally]
    push!(test_cases,
          Dict(:ode => ode,
               :funcs => funcs_to_test,
               :correct => Dict(funcs_to_test .=> correct)))

    #--------------------------------------------------------------------------

    ode = @ODEmodel(x1'(t)=-a1 * x1(t) + a21 * x2(t),
                    x2'(t)=-a2 * x2(t) - a21 * x2(t),
                    y1(t)=x1(t))
    funcs_to_test = [a1, a2, a21]
    correct = [:locally, :nonidentifiable, :nonidentifiable]
    push!(test_cases,
          Dict(:ode => ode,
               :funcs => funcs_to_test,
               :correct => Dict(funcs_to_test .=> correct)))

    #--------------------------------------------------------------------------

    ode = @ODEmodel(x1'(t)=-a1 * x1(t) + a21 * x2(t),
                    x2'(t)=-a2 * x2(t) - a21 * x2(t) + u(t),
                    y1(t)=x1(t))
    funcs_to_test = [a1, a2, a21, a2 + a1, a2 * (a1 - a21)]
    correct = [:locally, :locally, :globally, :globally, :globally]
    push!(test_cases,
          Dict(:ode => ode,
               :funcs => funcs_to_test,
               :correct => Dict(funcs_to_test .=> correct)))

    #--------------------------------------------------------------------------

    for case in test_cases
        result = assess_identifiability(case[:ode], case[:funcs])
        @test result == case[:correct]
    end
end
