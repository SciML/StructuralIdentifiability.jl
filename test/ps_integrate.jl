if GROUP == "All" || GROUP == "Core"
    @testset "Power series integration" begin
        T, t = Nemo.PowerSeriesRing(Nemo.QQ, 300, "t"; model = :capped_absolute)

        @test ps_integrate(zero(T)) == zero(T)

        @test ps_integrate(one(T)) == gen(T)

        @test ps_integrate(gen(T)) == gen(T)^2 * (1 // 2)

        @test ps_integrate(t^3 - t^4 * (1 // 8) + 5 * t^10) ==
              t^4 * (1 // 4) - t^5 * (1 // 40) + t^11 * (5 // 11)

        @test ps_integrate((1 + t)^1000) ==
              truncate((1 + t)^1001 * (1 // 1001) - 1 // 1001, 300)
    end
end
