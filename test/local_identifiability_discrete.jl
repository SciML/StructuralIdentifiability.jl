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

    cases = []

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

    for c in cases
        @test assess_local_identifiability(
            c[:dds];
            funcs_to_check = collect(keys(c[:res])),
            known_ic = c[:known],
        ) == c[:res]
    end
end
