@testset "Finding submodels" begin
    cases = [
        Dict(
            :ode => @ODEmodel(
                x1'(t) = x1(t)^2,
                x2'(t) = x1(t) * x2(t),
                y1(t) = x1(t),
                y2(t) = x2(t)
            ),
            :submodels => Set([(Set(["x1"]), Set(["y1"]), Set{String}())])
        ),
        Dict(
            :ode => @ODEmodel(
                x1'(t) = x1(t),
                x2'(t) = x2(t) + x1(t),
                y1(t) = x1(t),
                y2(t) = x1(t)^2,
                y3(t) = x2(t)
            ),
            :submodels => Set([(Set(["x1"]), Set(["y1", "y2"]), Set{String}())])
        ),
        Dict(
            :ode => @ODEmodel(
                x1'(t) = 0,
                x2'(t) = 0,
                x3'(t) = 0,
                y1(t) = x1(t) + x2(t),
                y2(t) = x2(t) + x3(t),
                y3(t) = x1(t) + x3(t)
            ),
            :submodels => Set([
                (Set(["x1", "x2"]), Set(["y1"]), Set{String}()),
                (Set(["x2", "x3"]), Set(["y2"]), Set{String}()),
                (Set(["x1", "x3"]), Set(["y3"]), Set{String}())
            ])
        ),
        Dict(
            :ode => @ODEmodel(
                x1'(t) = x1(t)^2,
                x2'(t) = x1(t) / x2(t),
                y1(t) = x1(t),
                y2(t) = x2(t)
            ),
            :submodels => Set([(Set(["x1"]), Set(["y1"]), Set{String}())])
        )
    ]

    for c in cases
        submodels = find_submodels(c[:ode])
        submodels = Set([(Set(map(var_to_str, ode.x_vars)), Set(map(var_to_str, ode.y_vars)), Set(map(var_to_str, ode.u_vars))) for ode in submodels])
        @test submodels == c[:submodels]
    end
end

