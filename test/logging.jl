if GROUP == "All" || GROUP == "Core"
    # An attempt to test logging
    using Logging

    @testset "Logging" begin
        ode = StructuralIdentifiability.@ODEmodel(x'(t) = x^42, y(t) = x)

        # Some logs
        @test_logs StructuralIdentifiability.assess_identifiability(ode)
        # No logs
        @test_logs min_level = Logging.Warn StructuralIdentifiability.assess_identifiability(
            ode,
            loglevel = Logging.Warn,
        )
        # Many logs
        @test_logs StructuralIdentifiability.assess_identifiability(
            ode,
            loglevel = Logging.Debug,
        )
    end
end
