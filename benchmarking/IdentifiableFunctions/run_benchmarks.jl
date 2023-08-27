# TODO: read benchmarks.jl, create and populate output directory from scratch

using ArgParse
using CpuId, Logging, Pkg, Printf
using Base.Threads
using Distributed
using Dates

using StructuralIdentifiability
using StructuralIdentifiability: _runtime_logger, ODE
using StructuralIdentifiability.ParamPunPam

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

include("benchmarks.jl")
include("utils.jl")

function parse_commandline()
    s = ArgParseSettings()
    #! format: off
    @add_arg_table s begin
        "--timeout"
            help = "Timeout, s."
            arg_type = Int
            default = 2500
        "--skip"
            help = "Systems to skip."
            arg_type = Vector{String}
            default = ["NFkB"]
        "--regen"
            help = "Re-generate the folder with benchmarks from scratch."
            arg_type = Bool
            default = true
        "--keywords"
            help = """
            Keyword arguments to `find_identifiable_functions`. 
            Semicolon-separated list of named tuples."""
            default = String
            default = "(strategy=(:gb, ),); (strategy=(:gb, ),with_states=true); (strategy=(:normalforms, 2),); (strategy=(:normalforms, 2),with_states=true); (strategy=(:normalforms, 3),); (strategy=(:normalforms, 3),with_states=true); (strategy=(:hybrid, ),); (strategy=(:hybrid, ),with_states=true)"
    end
    #! format: on

    parse_args(s)
end

function populate_benchmarks(args, kwargs)
    regen = args["regen"]
    !regen && return true
    @info "Re-generating the benchmarks folder"
    try
        if isdir((@__DIR__) * "/systems/")
            rm((@__DIR__) * "/systems/", recursive = true, force = true)
        end
    catch err
        @info "Something went wrong when deleting the benchmarks folder"
        showerror(stdout, err)
        println(stdout)
    end
    for bmark in benchmarks
        name = bmark[:name]
        system = bmark[:ode]
        @info "Generating $name"
        mkpath((@__DIR__) * "/systems/$name/")
        fd = open((@__DIR__) * "/systems/$name/$name.jl", "w")
        println(fd, "# $name")
        println(fd, "#! format: off")
        println(fd, "using StructuralIdentifiability")
        println(fd, "")
        ode = join(map(s -> "\t" * s, split(repr(system), "\n")[1:(end - 1)]), ",\n")
        println(fd, "system = @ODEmodel(\n$ode\n)")
        close(fd)
    end
    true
end

function run_benchmarks(args, kwargs)
    to_skip = args["skip"]
    timeout = args["timeout"]
    dirnames = first(walkdir((@__DIR__) * "/systems/"))[2]
    to_run_names = setdiff(dirnames, to_skip)
    to_run_indices = collect(1:length(to_run_names))

    nworkers = 4

    @info """
    Running benchmarks.
    Number of benchmarks: $(length(to_run_indices))
    Workers: $(nworkers)
    Timeout: $timeout seconds
    Keywords for `find_identifiable_functions`:
        $kwargs"""
    @info """
    Benchmark systems:
    $to_run_names"""

    seconds_passed(from_t) = round((time_ns() - from_t) / 1e9, digits = 2)

    queue = [(kw, idx) for kw in kwargs for idx in to_run_indices]
    procs = []
    log_fd = []
    keywords = []
    exited = []
    running = 0

    while true
        if !isempty(queue) && running < nworkers
            (kw, idx) = pop!(queue)
            name = to_run_names[idx]
            id = keywords_to_id(kw)
            @info "Running $name / $id"
            logs = open((@__DIR__) * "/systems/$name/logs_$id", "w")
            cmd =
                Cmd(["julia", (@__DIR__) * "/run_single_benchmark.jl", "$name", "$kw"])
            cmd = Cmd(cmd, ignorestatus = true, detach = false, env = copy(ENV))
            proc = run(pipeline(cmd, stdout = logs, stderr = logs), wait = false)
            push!(log_fd, logs)
            push!(keywords, kw)
            push!(
                procs,
                (index = idx, name = name, proc = proc, start_time = time_ns()),
            )
            running += 1
        end

        sleep(1.0)
        i = 1
        for i in 1:length(procs)
            i in exited && continue
            proc = procs[i]
            if process_exited(proc.proc)
                running -= 1
                push!(exited, i)
                close(log_fd[i])
                kw = keywords[i]
                start_time = proc.start_time
                @info "Yielded $(proc.name) / $(kw) after $(seconds_passed(start_time)) seconds"
            end
            if process_running(proc.proc)
                start_time = proc.start_time
                if seconds_passed(start_time) > timeout
                    kill(proc.proc)
                    close(log_fd[i])
                    kw = keywords[i]
                    running -= 1
                    push!(exited, i)
                    @warn "Timed-out $(proc.name) / $(kw) after $(seconds_passed(start_time)) seconds"
                end
            end
        end
        if length(exited) == length(to_run_names) * length(kwargs)
            @info "Exited $exited"
            @info "All benchmarks finished"
            break
        end
    end

    to_run_names
end

function collect_timings(args, kwargs, names; content = :compare)
    resulting_md = ""

    resulting_md *= """
    ## Simplification scores

    *Smaller is better.*

    $(now())

    Timeout: $(args["timeout"]) s
    """

    names = sort(names)
    runtimes = Dict()
    for name in names
        println("==== Reading $name")
        runtimes[name] = Dict()
        for kw in kwargs
            timings = nothing
            id = keywords_to_id(kw)
            runtimes[name][id] = Dict()
            try
                println("==== Opening /systems/$name/timings_$id")
                timings = open((@__DIR__) * "/systems/$name/timings_$id", "r")
            catch e
                @warn "Cannot collect timings for $name / $id"
                continue
            end
            lines = readlines(timings)
            if isempty(lines)
                @warn "Cannot collect timings for $name / $id"
                continue
            end
            @assert lines[1] == name
            for line in lines[2:end]
                k, v = split(line, ", ")
                runtimes[name][id][Symbol(k)] = parse(Float64, v)
            end
            close(timings)
        end
    end

    if content === :compare
        ids = map(keywords_to_id, kwargs)
        resulting_md *= "|Model|" * join(map(String ∘ Symbol, ids), "|") * "|\n"
        resulting_md *= "|-----|" * join(["---" for _ in ids], "|") * "|\n"
        for name in names
            times = runtimes[name]
            resulting_md *= "|$name|"
            for c in ids
                if isempty(times[c])
                    resulting_md *= " - " * "|"
                else
                    resulting_md *= @sprintf("%.2f", times[c][:id_total]) * "|"
                end
            end
            resulting_md *= "\n"
        end
    elseif length(content) == 2
        @assert content[1] === :compare
        feature = content[2]
        ids = map(keywords_to_id, kwargs)
        resulting_md *= "|Model|" * join(map(s -> String(Symbol(s)), ids), "|") * "|\n"
        resulting_md *= "|-----|" * join(["---" for _ in ids], "|") * "|\n"
        for name in names
            times = runtimes[name]
            resulting_md *= "|$name|"
            for c in ids
                if isempty(times[c])
                    resulting_md *= " - " * "|"
                else
                    # resulting_md *= @sprintf("%.2f", times[c][feature]) * "|"
                    resulting_md *= repr(round(Int, times[c][feature])) * "|"
                end
            end
            resulting_md *= "\n"
        end
    else
        kw = first(kwargs)
        id = keywords_to_id(kw)
        resulting_md *= "\nKeywords:\n$kw\n"
        resulting_md *=
            "|Model|" *
            join(map(c -> HUMAN_READABLE_CATEGORIES[c], ALL_CATEGORIES), "|") *
            "|\n"
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
    end

    resulting_md *= "\n*Benchmarking environment:*\n\n"
    resulting_md *= "* Total RAM (GiB): $(div(Sys.total_memory(), 2^30))\n"
    resulting_md *= "* Processor: $(cpubrand())\n"
    resulting_md *= "* Julia version: $(VERSION)\n\n"
    resulting_md *= "Versions of the dependencies:\n\n"

    deps = Pkg.dependencies()
    stid_info = deps[findfirst(x -> x.name == "StructuralIdentifiability", deps)]
    for (s, uid) in stid_info.dependencies
        if deps[uid].version !== nothing
            resulting_md *= "* $s : $(deps[uid].version)\n"
        end
    end

    open((@__DIR__) * "/benchmark_result.md", "w") do io
        write(io, resulting_md)
    end
end

function main()
    args = parse_commandline()
    kwargs = parse_keywords(args["keywords"])
    @info "Command-line args:"
    for (arg, val) in args
        @info "$arg  =>  $val"
    end
    @info "Keywords for `find_identifiable_functions`"
    @info kwargs
    flag = populate_benchmarks(args, kwargs)
    systems = run_benchmarks(args, kwargs)
    collect_timings(args, kwargs, systems, content = :compare)
end

main()
