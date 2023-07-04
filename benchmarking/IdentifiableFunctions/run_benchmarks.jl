using CpuId
using Logging
using Pkg
using Printf

using StructuralIdentifiability
using StructuralIdentifiability: _runtime_logger, ODE
using StructuralIdentifiability.ParamPunPam

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

include((@__DIR__) * "/benchmarks.jl")

runtimes = Dict()
ID_TIME_CATEGORIES = [
    :id_io_time,
    :id_global_time,
    :id_extract_funcs_time,
    :id_ideal_time,
    :id_filter_time,
    :id_inclusion_check,
    :id_groebner_time,
    :id_total,
]
GB_TIME_CATEGORIES = [
    :gb_discover_shape,
    :gb_discover_degrees,
    :gb_interpolate,
    :gb_recover_coeffs,
    :gb_npoints,
]
ALL_CATEGORIES = union(ID_TIME_CATEGORIES, GB_TIME_CATEGORIES)
NUM_RUNS = 1
NAME = ARGS[1]

# compile
find_identifiable_functions(benchmarks[1][:ode])

for bmark in benchmarks
    name = bmark[:name]
    if !(name == NAME)
        continue
    end
    if (:skip in keys(bmark)) && bmark[:skip]
        @info "Skipping $name"
        continue
    end
    @info "Processing $name"
    runtimes[bmark[:name]] = Dict(c => 0.0 for c in ALL_CATEGORIES)
    for _ in 1:NUM_RUNS
        funcs = find_identifiable_functions(bmark[:ode])
        @warn "===============\n$name\n$funcs\n==============="
        for cat in ID_TIME_CATEGORIES
            runtimes[name][cat] += _runtime_logger[cat]
        end
        for cat in GB_TIME_CATEGORIES
            runtimes[name][cat] += ParamPunPam._runtime_logger[cat]
        end
    end
    for k in keys(runtimes[name])
        runtimes[name][k] = runtimes[name][k] / NUM_RUNS
    end
    @warn "===============\n$name\n$runtimes\n==============="
end

resulting_md = ""

resulting_md *= "|Model|" * join(ALL_CATEGORIES, "|") * "|\n"
resulting_md *= "|-----|" * join(["---" for _ in ALL_CATEGORIES], "|") * "|\n"
for (name, times) in runtimes
    global resulting_md *= "|$name|"
    for c in ALL_CATEGORIES
        resulting_md *= @sprintf("%.2f", times[c]) * "|"
    end
    resulting_md *= "\n"
end

resulting_md *= "\n*Benchmarking environment:*\n\n"
resulting_md *= "* Total RAM (GiB): $(div(Sys.total_memory(), 2^30))\n"
resulting_md *= "* Processor: $(cpubrand())\n"
resulting_md *= "* Julia version: $(VERSION)\n\n"
resulting_md *= "Versions of the dependencies:\n\n"

deps = Pkg.dependencies()
stid_info = deps[findfirst(x -> x.name == "StructuralIdentifiability", deps)]
for (s, uid) in stid_info.dependencies
    if deps[uid].version != nothing
        global resulting_md *= "* $s : $(deps[uid].version)\n"
    end
end

open((@__DIR__) * "/benchmark_result.md", "w") do io
    write(io, resulting_md)
end
