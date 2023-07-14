using CpuId, Logging, Pkg, Printf
using Base.Threads
using Distributed
using Dates

using StructuralIdentifiability
using StructuralIdentifiability: _runtime_logger, ODE
using StructuralIdentifiability.ParamPunPam

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

ID_TIME_CATEGORIES = [
    :id_io_time,
    :id_global_time,
    :id_ideal_time,
    :id_inclusion_check,
    :id_groebner_time,
    :id_total,
]
ALL_CATEGORIES = ID_TIME_CATEGORIES

HUMAN_READABLE = Dict(
    :id_io_time => "io",
    :id_global_time => "global id.",
    :id_ideal_time => "gen. ideal",
    :id_inclusion_check => "inclusion",
    :id_groebner_time => "ParamPunPam.jl",
    :id_total => "total",
)

TO_SKIP = []
TO_RUN = []
TIMEOUT = 18_000

function run_benchmarks()
    dirnames = first(walkdir((@__DIR__) * "/systems/"))[2]
    to_run = isempty(TO_RUN) ? dirnames : intersect(dirnames, TO_RUN)
    to_run = setdiff(to_run, TO_SKIP)
    to_run_indices = collect(1:length(to_run))
    to_run_names = to_run
    @info """
    Running benchmarks.
    Number of benchmarks: $(length(to_run_indices))
    Workers: $(length(workers()))
    Timeout: $TIMEOUT sec."""
    @info """
    Benchmark systems:
    $to_run_names"""
    procs = []
    log_fd = []

    start_time = time_ns()
    for idx in to_run_indices
        name = to_run_names[idx]
        logs = open((@__DIR__) * "/systems/$name/logs", "w")
        cmd = Cmd(["julia", (@__DIR__) * "/run_single_benchmark.jl", "$name"])
        cmd = Cmd(cmd, ignorestatus = true, detach = false, env = copy(ENV))
        proc = run(pipeline(cmd, stdout = logs, stderr = logs), wait = false)
        push!(log_fd, logs)
        push!(procs, (index = idx, name = name, proc = proc))
    end

    ns_to_sec(start_time) = round((time_ns() - start_time) / 1e9, digits = 2)

    exited = []
    @label Wait
    sleep(1.0)
    i = 1
    for i in 1:length(procs)
        i in exited && continue
        proc = procs[i]
        if process_exited(proc.proc)
            push!(exited, i)
            close(log_fd[i])
            @info "Benchmark $(proc.name) yielded in $(ns_to_sec(start_time)) seconds"
        end
    end
    if ns_to_sec(start_time) > TIMEOUT
        @warn "Timed out after $(ns_to_sec(start_time)) seconds"
        for i in 1:length(procs)
            i in exited && continue
            kill(procs[i].proc)
            close(log_fd[i])
        end
    elseif length(exited) == length(procs)
        @info "All benchmarks finished"
    else
        @goto Wait
    end

    to_run_names
end

function collect_timings(names)
    resulting_md = ""

    resulting_md *= """
    ## Benchmark results

    $(now())


    """

    names = sort(names)
    runtimes = Dict()
    for name in names
        timings = nothing
        runtimes[name] = Dict()
        try
            timings = open((@__DIR__) * "/systems/$name/timings", "r")
        catch e
            @warn "Cannot collect timings for $name"
            continue
        end
        lines = readlines(timings)
        if isempty(lines)
            @warn "Cannot collect timings for $name"
            continue
        end
        @assert lines[1] == name
        for line in lines[2:end]
            k, v = split(line, ", ")
            runtimes[name][Symbol(k)] = parse(Float64, v)
        end
        close(timings)
    end

    resulting_md *=
        "|Model|" * join(map(c -> HUMAN_READABLE[c], ALL_CATEGORIES), "|") * "|\n"
    resulting_md *= "|-----|" * join(["---" for _ in ALL_CATEGORIES], "|") * "|\n"
    for name in names
        times = runtimes[name]
        resulting_md *= "|$name|"
        for c in ALL_CATEGORIES
            if isempty(times)
                resulting_md *= " - " * "|"
            else
                resulting_md *= @sprintf("%.2f", times[c]) * "|"
            end
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
            resulting_md *= "* $s : $(deps[uid].version)\n"
        end
    end

    open((@__DIR__) * "/benchmark_result.md", "w") do io
        write(io, resulting_md)
    end
end

names = run_benchmarks()
collect_timings(names)
