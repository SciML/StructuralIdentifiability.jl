# Adapter from https://github.com/sumiya11/Groebner.jl/blob/master/benchmark/CI-scripts/runtests.jl
# GPL license

# This script is not to be run directly. See runtests.jl.

stopwatch = time_ns()

# TTFX
t1 = @timed using StructuralIdentifiability

include("benchmarks.jl")

suite = []

push!(suite, (problem_name="using StructuralIdentifiability", type=:time, result=[t1.time]))

# Allocations
# Coming soon

# Runtime
import Primes
import Nemo

function nemo_make_prime_finite_field(p)
    if p < typemax(UInt)
        Nemo.fpField(convert(UInt, p), false)
    else
        Nemo.FpField(Nemo.ZZRingElem(p), false)
    end
end

function perform_operation(ode, fun=StructuralIdentifiability.assess_identifiability, trials=7; kws...)
    @asser fun in [
        StructuralIdentifiability.assess_identifiability,
        StructuralIdentifiability.assess_local_identifiability,
        StructuralIdentifiability.find_identifiable_functions,
        StructuralIdentifiability.reparametrize_global,
    ]
    times = []
    for _ in 1:trials
        GC.gc()
        time = @elapsed fun(ode; kws...)
        push!(times, time)
    end
    times
end

fun_name = "assess_identifiability"
perform_operation(benchmarks[:SIWR][:ode])
push!(
    suite,
    (
        problem_name=benchmarks[:SIWR][:name] * " " * fun_name,
        type=:time,
        result=perform_operation(benchmarks[:SIWR][:ode]),
    )
)

stopwatch = time_ns() - stopwatch
push!(suite, (problem_name="total", type=:time, result=[stopwatch / 1e9]))
