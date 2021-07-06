using CpuId
using Logging
using Pkg
using Printf

using StructuralIdentifiability
using StructuralIdentifiability: _runtime_logger, ODE

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

include("benchmarks.jl")

runtimes = Dict()
TIME_CATEGORIES = [:loc_time, :glob_time, :ioeq_time, :wrnsk_time, :rank_time, :check_time, :total]
NUM_RUNS = 5

for bmark in benchmarks
    name = bmark[:name]
    if (:skip in keys(bmark)) && bmark[:skip]
        @info "Skipping $name"
        continue
    end
    @info "Processing $name"
    runtimes[bmark[:name]] = Dict(c => 0. for c in TIME_CATEGORIES)
    for _ in 1:NUM_RUNS
        runtimes[name][:total] += @elapsed assess_identifiability(bmark[:ode])
        for cat in TIME_CATEGORIES[1:end - 1]
            runtimes[name][cat] += _runtime_logger[cat]
        end
    end
    for k in keys(runtimes[name])
        runtimes[name][k] = runtimes[name][k] / NUM_RUNS
    end
end

resulting_md = ""

resulting_md *= "|Model|" * join(TIME_CATEGORIES, "|") * "|\n"
resulting_md *= "|-----|" * join(["---" for _ in TIME_CATEGORIES], "|") * "|\n"
for (name, times) in runtimes
    global resulting_md *= "|$name|"
    for c in TIME_CATEGORIES
        resulting_md *= @sprintf("%.2f", times[c]) * "|"
    end
    resulting_md *= "\n"
end

resulting_md *= "\n*Benchmarking environment:*\n\n"
resulting_md *= "* Total RAM (Mb): $(Sys.total_memory() / 2^20)\n"
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


open("benchmark_result.md", "w") do io
    write(io, resulting_md)
end
