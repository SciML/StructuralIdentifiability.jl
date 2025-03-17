using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
Pkg.status()

Pkg.add(url = "https://github.com/SciML/StructuralIdentifiability.jl")

id = ARGS[1]

include("../common.jl")
include("../run_benchmarks.jl")

dump_results((@__DIR__) * "/results_$id", "master")
