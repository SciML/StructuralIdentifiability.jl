@testset "Finding submodels" begin
    cases = [
        Dict(
            :ode => @ODEmodel(
                x1'(t) = x1(t)^2,
                x2'(t) = x1(t) * x2(t),
                y1(t) = x1(t),
                y2(t) = x2(t)
            ),
            :submodels => Set([(["x1(t)"], ["y1(t)"], Vector{String}())]),
        ),
        Dict(
            :ode => @ODEmodel(
                x1'(t) = x1(t),
                x2'(t) = x2(t) + x1(t),
                y1(t) = x1(t),
                y2(t) = x1(t)^2,
                y3(t) = x2(t)
            ),
            :submodels => Set([(["x1(t)"], ["y1(t)", "y2(t)"], Vector{String}())]),
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
                (["x1(t)", "x2(t)"], ["y1(t)"], Vector{String}()),
                (["x2(t)", "x3(t)"], ["y2(t)"], Vector{String}()),
                (["x1(t)", "x3(t)"], ["y3(t)"], Vector{String}()),
            ]),
        ),
        Dict(
            :ode => @ODEmodel(
                x1'(t) = x1(t)^2,
                x2'(t) = x1(t) / x2(t),
                y1(t) = x1(t),
                y2(t) = x2(t)
            ),
            :submodels => Set([(["x1(t)"], ["y1(t)"], Vector{String}())]),
        ),
        Dict(
            :ode => @ODEmodel(
                x1'(t) = x1(t) + a(t) * x2(t)^2,
                x2'(t) = x2(t) - b(t) * x1(t) * x2(t),
                x3'(t) = x3(t) + x1(t) + x2(t),
                y1(t) = d(t) * x1(t),
                y2(t) = x3(t)
            ),
            :submodels =>
                Set([(["x1(t)", "x2(t)"], ["y1(t)"], ["a(t)", "b(t)", "d(t)"])]),
        ),
        Dict(
            :ode => @ODEmodel(
                x1'(t) = x1(t),
                x2'(t) = x1(t) - c * a(t) * x2(t) * x2(t),
                x3'(t) = x3(t) - x1(t) * x2(t) * x3(t),
                y1(t) = x1(t),
                y2(t) = x1(t) + d(t) * x2(t),
                y3(t) = x3(t)
            ),
            :submodels => Set([
                (["x1(t)", "x2(t)"], ["y1(t)", "y2(t)"], ["a(t)", "d(t)"]),
                (["x1(t)"], ["y1(t)"], Vector{String}()),
            ]),
        ),
        Dict(
            :ode => @ODEmodel(
                o1'(t) = r1 * o1(t),
                x0'(t) = o1(t) - d * x0(t),
                x1'(t) = x0(t) - x1(t),
                x2'(t) = x1(t) - x2(t),
                x3'(t) = x2(t) - x3(t),
                x4'(t) = x3(t) - x4(t),
                x5'(t) = x4(t) - x5(t),
                o12'(t) = r2 * o12(t),
                x02'(t) = o12(t) + x2(t) + x3(t) + x22(t) + x32(t) + x02(t),
                x12'(t) = x2(t) + x3(t) + x22(t) + x32(t) + x02(t) - x12(t),
                x22'(t) = x12(t) - x22(t),
                x32'(t) = x22(t) - x32(t),
                x42'(t) = x32(t) - x42(t),
                x52'(t) = x42(t) - x52(t),
                y0(t) = x0(t) + x02(t),
                y1(t) = x1(t) + x12(t),
                y2(t) = x2(t) + x22(t),
                y3(t) = x3(t) + x32(t),
                y4(t) = x4(t) + x42(t) + x5(t) + x52(t),
                y5(t) =
                    x1(t) + x2(t) + x3(t) + x4(t) + x12(t) + x22(t) + x32(t) + x42(t),
                y6(t) = x5(t) + x52(t)
            ),
            :submodels => Set([
                (
                    [
                        "o1(t)",
                        "x0(t)",
                        "x1(t)",
                        "x2(t)",
                        "x3(t)",
                        "o12(t)",
                        "x02(t)",
                        "x12(t)",
                        "x22(t)",
                        "x32(t)",
                    ],
                    ["y0(t)", "y1(t)", "y2(t)", "y3(t)"],
                    Vector{String}(),
                ),
                (
                    [
                        "o1(t)",
                        "x0(t)",
                        "x1(t)",
                        "x2(t)",
                        "x3(t)",
                        "x4(t)",
                        "o12(t)",
                        "x02(t)",
                        "x12(t)",
                        "x22(t)",
                        "x32(t)",
                        "x42(t)",
                    ],
                    ["y0(t)", "y1(t)", "y2(t)", "y3(t)", "y5(t)"],
                    Vector{String}(),
                ),
            ]),
        ),
    ]

    for c in cases
        submodels = find_submodels(c[:ode])
        submodels = Set([
            (
                collect(map(var_to_str, ode.x_vars)),
                collect(map(var_to_str, ode.y_vars)),
                collect(map(var_to_str, ode.u_vars)),
            ) for ode in submodels
        ])
        @test submodels == c[:submodels]
    end
end
