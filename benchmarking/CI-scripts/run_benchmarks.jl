# Adapter from https://github.com/sumiya11/Groebner.jl/blob/master/benchmark/CI-scripts/runtests.jl
# GPL license

# This script is not to be run directly. See runtests.jl.

stopwatch = time_ns()

# TTFX
t1 = @timed using StructuralIdentifiability

suite = []

push!(
    suite,
    (problem_name = "using StructuralIdentifiability", type = :time, result = [t1.time]),
)

# Record the load time before @ODEmodel introduces model parameters such as t1.
include("../benchmarks.jl")

# Allocations
# Coming soon

# Runtime
import Primes
import Nemo
using Logging

function nemo_make_prime_finite_field(p)
    return if p < typemax(UInt)
        Nemo.fpField(convert(UInt, p), false)
    else
        Nemo.FpField(Nemo.ZZRingElem(p), false)
    end
end

function perform_operation(
        ode,
        fun = StructuralIdentifiability.assess_identifiability,
        trials = 3;
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
        time = @elapsed fun(ode; loglevel=Logging.Error, kws...)
        push!(times, time)
    end
    return times
end

fun_name = "assess_identifiability"
perform_operation(benchmarks[:LV][:ode])
for (m, trials) in [(:SIWR_orig, 3), (:LV, 30), (:Pharm, 3)]
    push!(
        suite,
        (
            problem_name = benchmarks[m][:name] * " " * fun_name,
            type = :time,
            result = perform_operation(
                benchmarks[m][:ode],
                StructuralIdentifiability.assess_identifiability,
                trials,
            ),
        ),
    )
end

fun_name = "assess_local_identifiability"
perform_operation(
    benchmarks[:LV][:ode],
    StructuralIdentifiability.assess_local_identifiability,
)
for (m, trials) in [
        (:SIWR_orig, 30), (:Pharm, 30), (:MAPK_5out, 3), (:MAPK_5out_bis, 3), (:Goodwin, 30),
    ]
    push!(
        suite,
        (
            problem_name = benchmarks[m][:name] * " " * fun_name,
            type = :time,
            result = perform_operation(
                benchmarks[m][:ode],
                StructuralIdentifiability.assess_local_identifiability,
                trials,
            ),
        ),
    )
end

fun_name = "find_identifiable_functions"
perform_operation(
    benchmarks[:LV][:ode],
    StructuralIdentifiability.find_identifiable_functions,
)
for (m, trials) in [
        (:Goodwin, 30), (:SEAIJRC, 3), (:St, 3), (:LLW, 30), (:Bilirubin, 30),
    ]
    push!(
        suite,
        (
            problem_name = benchmarks[m][:name] * " " * fun_name,
            type = :time,
            result = perform_operation(
                benchmarks[m][:ode],
                StructuralIdentifiability.find_identifiable_functions,
                trials,
            ),
        ),
    )
end

fun_name = "reparametrize_global"
perform_operation(
    benchmarks[:LV][:ode],
    StructuralIdentifiability.reparametrize_global,
)
for (m, trials) in [(:Goodwin, 30), (:SEAIJRC, 3), (:SEUIR, 30)]
    push!(
        suite,
        (
            problem_name = benchmarks[m][:name] * " " * fun_name,
            type = :time,
            result = perform_operation(
                benchmarks[m][:ode],
                StructuralIdentifiability.reparametrize_global,
                trials,
            ),
        ),
    )
end

stopwatch = time_ns() - stopwatch
push!(suite, (problem_name = "total", type = :time, result = [stopwatch / 1.0e9]))
