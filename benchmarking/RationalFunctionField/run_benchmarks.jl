using Dates
using Printf
using BenchmarkTools

using Nemo

using Logging

using StructuralIdentifiability: simplified_generating_set, RationalFunctionField, restart_logging

cases_simplification = []

cases_containment = []

include("from_nouvelle_vague.jl")
include("from_io_equations.jl")

resulting_md = """
## Timings for simplification

|Case description|weak|standard|
|----------------|----|--------|
"""

restart_logging(loglevel = Logging.Debug)

function run_case(case)
    @warn "Starting " * case[:description]
    g = case[:gens]
    description = case[:description]
    runtimes = []
    for level in (:weak, :standard)
        runtime = @belapsed simplified_generating_set(RationalFunctionField($g), simplify = $level)
        push!(runtimes, runtime)
    end
    return runtimes
end

for case in cases_simplification
    runtimes = run_case(case)
    global resulting_md *= "|" * case[:description] * "|" * join(map(t -> @sprintf("%.2f", t), runtimes), "|") * "|\n"
end

open("benchmark_results_$(today()).md", "w") do io
    write(io, resulting_md)
end
