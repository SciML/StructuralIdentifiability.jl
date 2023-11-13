@testset "Assessing identifiability" begin
    test_cases = []

    # 2-compartiment model

    ode = @ODEmodel(
        x0'(t) = -(a01 + a21) * x0(t) + a12 * x1(t),
        x1'(t) = a21 * x0(t) - a12 * x1(t),
        y(t) = x0(t)
    )
    funcs_to_test =
        [a01, a21, a12, a01 * a12, a01 + a12 + a21, (a01 + a12 + a21) // (a01 * a12)]
    correct = [
        :nonidentifiable,
        :nonidentifiable,
        :nonidentifiable,
        :globally,
        :globally,
        :globally,
    ]
    push!(
        test_cases,
        Dict(
            :ode => ode,
            :funcs => funcs_to_test,
            :correct => Dict(funcs_to_test .=> correct),
        ),
    )

    #--------------------------------------------------------------------------
    # No parameters no worry

    ode = @ODEmodel(x1'(t) = x1, x2'(t) = x2, y(t) = x1 + x2(t))
    funcs_to_test = [x1, x2, x1 + x2]
    correct = [:nonidentifiable, :nonidentifiable, :globally]
    push!(
        test_cases,
        Dict(
            :ode => ode,
            :funcs => funcs_to_test,
            :correct => Dict(funcs_to_test .=> correct),
        ),
    )

    # Also test when `funcs_to_test` is empty!
    funcs_to_test = Vector{typeof(x1)}()
    correct = Dict(x1 => :nonidentifiable, x2 => :nonidentifiable)
    push!(test_cases, Dict(:ode => ode, :funcs => funcs_to_test, :correct => correct))

    #--------------------------------------------------------------------------

    ode = @ODEmodel(
        x0'(t) = a * x0(t) - b * x0(t) * x1(t) + u(t),
        x1'(t) = c * x1(t) + d * x0(t) * x1(t),
        y(t) = x0(t)
    )
    funcs_to_test = [a, b, c, d, b * x1, x0, x1]
    correct = [
        :globally,
        :nonidentifiable,
        :globally,
        :globally,
        :globally,
        :globally,
        :nonidentifiable,
    ]
    push!(
        test_cases,
        Dict(
            :ode => ode,
            :funcs => funcs_to_test,
            :correct => Dict(funcs_to_test .=> correct),
        ),
    )

    #--------------------------------------------------------------------------

    ode = @ODEmodel(
        S'(t) = mu - bi * S(t) * I(t) - bw * S(t) * W(t) - mu * S(t) + a * R(t),
        I'(t) = bw * S(t) * W(t) + bi * S(t) * I(t) - (gam + mu) * I(t),
        W'(t) = xi * (I(t) - W(t)),
        R'(t) = gam * I(t) - (mu + a) * R(t),
        y(t) = k * I(t)
    )
    funcs_to_test = [mu, bi, bw, a, xi, gam, mu, gam + mu, k]
    correct = [:globally for _ in funcs_to_test]
    push!(
        test_cases,
        Dict(
            :ode => ode,
            :funcs => funcs_to_test,
            :correct => Dict(funcs_to_test .=> correct),
        ),
    )

    #--------------------------------------------------------------------------

    ode = @ODEmodel(
        x1'(t) = -b * x1(t) + 1 / (c + x4(t)),
        x2'(t) = alpha * x1(t) - beta * x2(t),
        x3'(t) = gama * x2(t) - delta * x3(t),
        x4'(t) = sigma * x4(t) * (gama * x2(t) - delta * x3(t)) / x3(t),
        y(t) = x1(t)
    )
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
    push!(
        test_cases,
        Dict(
            :ode => ode,
            :funcs => funcs_to_test,
            :correct => Dict(funcs_to_test .=> correct),
        ),
    )

    #--------------------------------------------------------------------------

    ode = @ODEmodel(
        x1'(t) = (1 + x1(t)^2) // 2,
        x2'(t) = (1 - x1(t)^2) // (1 + x1(t)^2),
        y1(t) = 2 * x1(t) // (b * (1 + x1(t)^2)),
        y2(t) = x2(t)
    )
    funcs_to_test = [b, x1, x2]
    correct = [:globally, :globally, :globally]
    push!(
        test_cases,
        Dict(
            :ode => ode,
            :funcs => funcs_to_test,
            :correct => Dict(funcs_to_test .=> correct),
        ),
    )

    #--------------------------------------------------------------------------

    ode = @ODEmodel(
        x1'(t) = -a1 * x1(t) + a21 * x2(t),
        x2'(t) = -a2 * x2(t) - a21 * x2(t),
        y1(t) = x1(t)
    )
    funcs_to_test = [a1, a2, a21]
    correct = [:locally, :nonidentifiable, :nonidentifiable]
    push!(
        test_cases,
        Dict(
            :ode => ode,
            :funcs => funcs_to_test,
            :correct => Dict(funcs_to_test .=> correct),
        ),
    )

    #--------------------------------------------------------------------------

    ode = @ODEmodel(
        x1'(t) = -a1 * x1(t) + a21 * x2(t),
        x2'(t) = -a2 * x2(t) - a21 * x2(t) + u(t),
        y1(t) = x1(t)
    )
    funcs_to_test = [a1, a2, a21, a2 + a1, a2 * (a1 - a21)]
    correct = [:locally, :locally, :globally, :globally, :globally]
    push!(
        test_cases,
        Dict(
            :ode => ode,
            :funcs => funcs_to_test,
            :correct => Dict(funcs_to_test .=> correct),
        ),
    )

    #-------------------------------------------------------------------------- 

    ode = @ODEmodel(
        x1'(t) = -(a21 + a31 + a01) * x1(t) + a12 * x2(t) + a13 * x3(t) + u(t),
        x2'(t) = a21 * x1(t) - a12 * x2(t),
        x3'(t) = a31 * x1(t) - a13 * x3(t),
        y(t) = x1(t)
    )
    funcs_to_test = [x1, x2, x3, a31 + a21, (a21 - a31) // (a12 - a13)]
    correct = [:globally, :locally, :locally, :globally, :globally]
    push!(
        test_cases,
        Dict(
            :ode => ode,
            :funcs => funcs_to_test,
            :correct => Dict(funcs_to_test .=> correct),
        ),
    )

    #-------------------------------------------------------------------------- 

    ode = @ODEmodel(
        S'(t) = -b * S(t) * In(t) / N,
        E'(t) = b * S(t) * In(t) / N - nu * E(t),
        In'(t) = nu * E(t) - a * In(t),
        y1(t) = In(t),
        y2(t) = N
    )
    funcs_to_test = [b, N, In, a, nu, S, E, a * nu, a + nu]
    correct = [
        :globally,
        :globally,
        :globally,
        :locally,
        :locally,
        :locally,
        :locally,
        :globally,
        :globally,
    ]
    push!(
        test_cases,
        Dict(
            :ode => ode,
            :funcs => funcs_to_test,
            :correct => Dict(funcs_to_test .=> correct),
        ),
    )

    #-------------------------------------------------------------------------- 

    for case in test_cases
        result = assess_identifiability(case[:ode], funcs_to_check = case[:funcs])
        @test result == case[:correct]
    end
end
