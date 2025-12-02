# Test that models in benchmarking/benchmarks.jl are well-formed.

benchmarking_dir = joinpath(dirname(dirname(pathof(StructuralIdentifiability))), "benchmarking")
include(joinpath(benchmarking_dir, "benchmarks.jl"))

@testset "Benchmark definitions" begin
    for (name, bench) in benchmarks
        @test haskey(bench, :name)
        @test haskey(bench, :ode)
        if haskey(bench, :cite)
            @test isfile(joinpath(benchmarking_dir, bench[:cite]))
        end
    end
end
