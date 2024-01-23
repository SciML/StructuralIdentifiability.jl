if GROUP == "All" || GROUP == "Core"
    @testset "Discrete local identifiability, internal function" begin
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

        dds = @DDSmodel(a(t + 1) = b(t) * a(t) + c, b(t + 1) = d * a(t), y(t) = b(t))

        push!(
            cases,
            Dict(
                :dds => dds,
                :res => OrderedDict(
                    b => true,
                    a => false,
                    c => false,
                    d => false,
                    d * a => true,
                    d * c => true,
                ),
                :known => :none,
            ),
        )

        push!(
            cases,
            Dict(
                :dds => dds,
                :res => OrderedDict(b => true, a => true, c => true, d => true),
                :known => [a],
            ),
        )

        push!(
            cases,
            Dict(
                :dds => dds,
                :res => OrderedDict(b => true, a => true, c => true, d => true),
                :known => :all,
            ),
        )

        # -------------------

        # Example 4 from https://doi.org/10.1016/j.automatica.2016.01.054
        dds = @DDSmodel(x(t + 1) = theta^3 * x(t), y(t) = x(t))

        push!(
            cases,
            OrderedDict(
                :dds => dds,
                :res => Dict(theta => true, x => true),
                :known => :none,
            ),
        )

        # -------------------

        for c in cases
            @test _assess_local_identifiability_discrete_aux(
                c[:dds],
                collect(keys(c[:res])),
                c[:known],
            ) == c[:res]
        end
    end
end
