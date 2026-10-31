include(joinpath(@__DIR__, "..", "shared", "test_setup.jl"))
include(joinpath(@__DIR__, "../..", "benchmarking", "benchmarks.jl"))

if GROUP == "All" || GROUP == "Core"
    @testset "Output saturation" begin
        ode = @ODEmodel(
            x'(t) = z(t),
            z'(t) = 0,
            y(t) = x(t) // a
        )
        saturated = saturate_outputs(
            ode,
            Dict(only(ode.y_vars) => 1),
        )
        @test length(saturated.y_vars) == 2


        # full nfkb model from benchmarks
        nfkb = benchmarks[:NFkB_full][:ode]

        nfkb_orders = propose_orders(nfkb)
        new_nfkb = saturate_outputs(nfkb, nfkb_orders)
        @test length(new_nfkb.y_vars) == 14
        nfkb_time = @elapsed nfkb_ident_results = assess_identifiability(nfkb)
        @test nfkb_time < 500
        @test count(x -> x == :globally, values(nfkb_ident_results)) == 37
        @test count(x -> x == :nonidentifiable, values(nfkb_ident_results)) == 7

        # TumorPilis model from the benchmarks
        tum_pil = benchmarks[:TumorPillis][:ode]

        tum_pil_orders = propose_orders(tum_pil)
        new_tum_pil = saturate_outputs(tum_pil, tum_pil_orders)
        @test length(new_tum_pil.y_vars) == 6
        tum_pil_time = @elapsed tum_pil_results = assess_identifiability(tum_pil)
        @test tum_pil_time < 500
        @test count(x -> x == :globally, values(tum_pil_results)) == 22
        @test count(x -> x == :nonidentifiable, values(tum_pil_results)) == 9

        weckel = benchmarks[:Weckel][:ode]
        weckel_orders = propose_orders(weckel)
        new_weckel = saturate_outputs(weckel, weckel_orders)
        @test length(new_weckel.y_vars) == 10
        weckel_time = @elapsed weckel_results = assess_identifiability(weckel)
        @test weckel_time < 500
        @test count(x -> x == :globally, values(weckel_results)) == 12
        @test count(x -> x == :nonidentifiable, values(weckel_results)) == 6
    end
end
