@testset "Discrete local identifiability, @DDSmodel interface" begin
    cases = []

    dds = @DDSmodel(a(t + 1) = (b + c) * a(t) + 1, y(t) = a(t))

    push!(
        cases,
        Dict(
            :dds => dds,
            :res => OrderedDict(a => true, b => false, c => false, b + c => true),
            :known => :none,
        ),
    )

    push!(
        cases,
        Dict(
            :dds => dds,
            :res => OrderedDict(a => true, b => false, c => false, b + c => true),
            :known => :all,
        ),
    )

    #---------------------

    dds = @DDSmodel(x1(t + 1) = a * x1(t) * x2(t), x2(t + 1) = b * x1(t), y(t) = x2(t))

    push!(
        cases,
        Dict(
            :dds => dds,
            :res =>
                OrderedDict(a => true, b => false, x1 => false, x2 => true, b * x1 => true),
            :known => :none,
        ),
    )

    push!(
        cases,
        Dict(
            :dds => dds,
            :res => OrderedDict(a => true, b => true, x1 => true, x2 => true),
            :known => :all,
        ),
    )

    #---------------------

    # Example 1 from https://doi.org/10.1016/j.automatica.2008.03.019
    dds = @DDSmodel(
        x1(t + 1) = theta1 * x1(t) + x2(t),
        x2(t + 1) = (1 - theta2) * x1(t) + x2(t)^2 + u(t),
        y(t) = x1(t)
    )

    push!(
        cases,
        Dict(
            :dds => dds,
            :res => OrderedDict(x1 => true, x2 => true, theta1 => true, theta2 => true),
            :known => :none,
        ),
    )

    # Example 2 from https://doi.org/10.1016/j.automatica.2008.03.019
    dds = @DDSmodel(
        x1(t + 1) = theta1 * x1(t)^2 + theta2 * x2(t) + u(t),
        x2(t + 1) = theta3 * x1(t),
        y(t) = x1(t)
    )

    push!(
        cases,
        Dict(
            :dds => dds,
            :res => OrderedDict(
                x1 => true,
                x2 => false,
                theta1 => true,
                theta2 => false,
                theta3 => false,
            ),
            :known => :none,
        ),
    )

    dds = @DDSmodel(
        x1(t + 1) = theta1 * x1(t)^2 + theta2 * x2(t) + u(t),
        x2(t + 1) = theta3 * x1(t),
        y(t) = x1(t),
        y2(t) = x2(t)
    )
    push!(
        cases,
        Dict(
            :dds => dds,
            :res => Dict(
                x1 => true,
                x2 => true,
                theta1 => true,
                theta2 => true,
                theta3 => true,
            ),
            :known => :none,
        ),
    )

    #---------------------

    for c in cases
        @test assess_local_identifiability(
            c[:dds];
            funcs_to_check = collect(keys(c[:res])),
            known_ic = c[:known],
        ) == c[:res]
    end
end
