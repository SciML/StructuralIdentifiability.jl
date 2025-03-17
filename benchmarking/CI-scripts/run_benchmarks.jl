# Adapter from https://github.com/sumiya11/Groebner.jl/blob/master/benchmark/CI-scripts/runtests.jl
# GPL license

# This script is not to be run directly. See runtests.jl.

stopwatch = time_ns()

# TTFX
t1 = @timed using StructuralIdentifiability

include("benchmarks.jl")

suite = []

push!(
    suite,
    (problem_name = "using StructuralIdentifiability", type = :time, result = [t1.time]),
)

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

function perform_operation(
    ode,
    fun = StructuralIdentifiability.assess_identifiability,
    trials = 5;
    kws...,
)
    @assert fun in [
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
for m in [:SIWR, :LV_simple, :Pharm, :MAPK6, :Goodwin]
    push!(
        suite,
        (
            problem_name = benchmarks[m][:name] * " " * fun_name,
            type = :time,
            result = perform_operation(benchmarks[m][:ode]),
        ),
    )
end

fun_name = "assess_local_identifiability"
perform_operation(benchmarks[:SIWR][:ode])
for m in [:SIWR, :Pharm, :MAPK5, :MAPK5bis, :Goodwin]
    push!(
        suite,
        (
            problem_name = benchmarks[m][:name] * " " * fun_name,
            type = :time,
            result = perform_operation(
                benchmarks[m][:ode],
                StructuralIdentifiability.assess_local_identifiability,
            ),
        ),
    )
end

fun_name = "find_identifiable_functions"
perform_operation(benchmarks[:SIWR][:ode])
for m in [:SIWR, :Goodwin, :SEAIJRC, :Sntg, :QY, :LLW, :Bilirubin]
    push!(
        suite,
        (
            problem_name = benchmarks[m][:name] * " " * fun_name,
            type = :time,
            result = perform_operation(
                benchmarks[m][:ode],
                StructuralIdentifiability.find_identifiable_functions,
            ),
        ),
    )
end

fun_name = "reparametrize_global"
perform_operation(benchmarks[:SIWR][:ode])
for m in [:Goodwin, :SEAIJRC, :Sntg, :LLW, :Bilirubin, :SEUIR]
    push!(
        suite,
        (
            problem_name = benchmarks[m][:name] * " " * fun_name,
            type = :time,
            result = perform_operation(
                benchmarks[m][:ode],
                StructuralIdentifiability.reparametrize_global,
            ),
        ),
    )
end

stopwatch = time_ns() - stopwatch
push!(suite, (problem_name = "total", type = :time, result = [stopwatch / 1e9]))
