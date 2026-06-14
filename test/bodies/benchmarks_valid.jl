include(joinpath(@__DIR__, "..", "shared", "test_setup.jl"))

@testset "Benchmarks are valid" verbose = true begin
    # https://github.com/pogudingleb/RationalFunctionFields.jl uses this
    include(joinpath(dirname(dirname(pathof(StructuralIdentifiability))), "benchmarking", "benchmarks.jl"))
end
