using ArgParse
using CpuId, Logging, Pkg, Printf
using Base.Threads
using Distributed
using Dates
using ProgressMeter

using StructuralIdentifiability
using StructuralIdentifiability: _runtime_logger, ODE
using StructuralIdentifiability.ParamPunPam

global_logger(Logging.SimpleLogger(stdout, Logging.Warn))
include("benchmarks.jl")
global_logger(Logging.SimpleLogger(stdout, Logging.Info))

const _progressbar_color = :light_green
const _progressbar_value_color = :light_green
progressbar_enabled() =
    Logging.Info <= Logging.min_enabled_level(current_logger()) < Logging.Warn

include("utils.jl")

function parse_commandline()
    s = ArgParseSettings()
    #! format: off
    @add_arg_table s begin
        "--timeout"
            help = "Timeout, s."
            arg_type = Int
            default = 300
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
        "--reparam"
            help = "Do model reparametrization."
            arg_type = Bool
            default = false
    end
    #! format: on

    parse_args(s)
end

function populate_benchmarks(args, kwargs)
    regen = args["regen"]
    !regen && return true
    @debug "Re-generating the benchmarks folder"
    try
        if isdir((@__DIR__) * "/systems/")
            rm((@__DIR__) * "/systems/", recursive = true, force = true)
        end
    catch err
        @info "Something went wrong when deleting the benchmarks folder"
        showerror(stdout, err)
        println(stdout)
    end
    prog = Progress(
        length(benchmarks),
        "Generating benchmarks",
        # spinner = true,
        dt = 0.1,
        enabled = progressbar_enabled(),
        color = _progressbar_color,
    )
    for bmark in benchmarks
        next!(prog) # , spinner = "⌜⌝⌟⌞")
        name = bmark[:name]
        system = bmark[:ode]
        @debug "Generating $name"
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
    finish!(prog)
    true
end

function make_tasks(reparam, problems, kwargs)
    if reparam
        [(problem_id = problem, cli_args = NamedTuple()) for problem in problems]
    else
        [(problem_id = problem, cli_args = kw) for kw in kwargs for problem in problems]
    end
end

function run_task(reparam) end

function run_benchmarks(args, kwargs)
    to_skip = args["skip"]
    timeout = args["timeout"]
    reparam = args["reparam"]
    dirnames = first(walkdir((@__DIR__) * "/systems/"))[2]
    to_run_names = setdiff(dirnames, to_skip)
    to_run_indices = collect(1:length(to_run_names))

    nworkers = 16

    if reparam
        @info """
        Model reparametrization."""
    else
        @info """
        Identifiable functions.

        Keywords for `find_identifiable_functions`:
        \t$(join(map(string, kwargs), "\n\t"))"""
    end

    @info """
    Running benchmarks.
    Number of benchmark systems: $(length(to_run_indices))
    Workers: $(nworkers)
    Timeout: $timeout seconds"""
    @info """
    Benchmark systems:
    $to_run_names"""

    seconds_passed(from_t) = round((time_ns() - from_t) / 1e9, digits = 2)

    queue = make_tasks(reparam, to_run_indices, kwargs)
    procs = []
    log_fd = []
    keywords = []
    exited = []
    errored = []
    running = 0

    generate_showvalues(procs) =
        () -> [(
            :Active,
            join(
                map(
                    proc -> string(proc.name) * " / " * string(proc.global_id),
                    filter(proc -> process_running(proc.proc), procs),
                ),
                ", ",
            ),
        )]

    prog = Progress(
        length(queue),
        "Running benchmarks",
        # spinner = true,
        dt = 0.3,
        enabled = progressbar_enabled(),
        color = _progressbar_color,
    )

    while true
        if !isempty(queue) && running < nworkers
            next!(
                prog,
                showvalues = generate_showvalues(procs),
                step = 0,
                valuecolor = _progressbar_value_color,
                # spinner = "⌜⌝⌟⌞",
            )
            task = pop!(queue)
            cli_args = task.cli_args
            problem_idx = to_run_indices[task.problem_id]
            problem_name = to_run_names[problem_idx]
            global_id = keywords_to_id(cli_args)
            @debug "Running $name / $id"
            logs = open((@__DIR__) * "/systems/$problem_name/logs_$global_id", "w")
            if reparam
                cmd = Cmd([
                    "julia",
                    (@__DIR__) * "/run_single_benchmark_reparam.jl",
                    "$problem_name",
                    # "$kw",
                ])
            else
                cmd = Cmd([
                    "julia",
                    (@__DIR__) * "/run_single_benchmark.jl",
                    "$problem_name",
                    "$kw",
                ])
            end
            cmd = Cmd(cmd, ignorestatus = true, detach = false, env = copy(ENV))
            proc = run(pipeline(cmd, stdout = logs, stderr = logs), wait = false)
            push!(log_fd, logs)
            push!(keywords, cli_args)
            push!(
                procs,
                (
                    index = problem_idx,
                    name = problem_name,
                    proc = proc,
                    start_time = time_ns(),
                    global_id = global_id,
                ),
            )
            running += 1
        end

        sleep(0.2)
        i = 1
        for i in 1:length(procs)
            i in exited && continue
            proc = procs[i]
            if process_exited(proc.proc)
                running -= 1
                push!(exited, i)
                if proc.proc.exitcode != 0
                    push!(errored, i)
                end
                close(log_fd[i])
                kw = keywords[i]
                start_time = proc.start_time
                next!(
                    prog,
                    showvalues = generate_showvalues(procs),
                    valuecolor = _progressbar_value_color,
                )
                @debug "Yielded $(proc.name) / $(kw) after $(seconds_passed(start_time)) seconds"
            end
            if process_running(proc.proc)
                start_time = proc.start_time
                if seconds_passed(start_time) > timeout
                    kill(proc.proc)
                    close(log_fd[i])
                    kw = keywords[i]
                    running -= 1
                    push!(exited, i)
                    next!(
                        prog,
                        showvalues = generate_showvalues(procs),
                        valuecolor = _progressbar_value_color,
                    )
                    @debug "Timed-out $(proc.name) / $(kw) after $(seconds_passed(start_time)) seconds"
                end
            end
        end
        if length(exited) == length(to_run_names) * 1 # * length(kwargs)
            @debug "Exited $exited"
            @debug "All benchmarks finished"
            break
        end
    end
    finish!(prog)

    if !isempty(errored)
        printstyled("(!) Maybe errored:\n", color = :red)
        for i in errored
            println("\t$(procs[i].name) / $(procs[i].global_id)")
        end
    end

    to_run_names
end

function collect_timings(args, kwargs, names; content = :compare)
    resulting_md = ""

    kwargs = [[]]

    resulting_md *= """
    ## Benchmark results

    Timestamp: $(now())

    Timeout: $(args["timeout"]) s

    **Timings in seconds.**

    """

    cannot_collect = []
    names = sort(names)
    runtimes = Dict()
    data = Dict()
    for name in names
        @debug "==== Reading $name"
        runtimes[name] = Dict()
        data[name] = Dict()
        for kw in kwargs
            timings = nothing
            id = keywords_to_id(kw)
            #####
            runtimes[name][id] = Dict()
            try
                @debug "==== Opening /systems/$name/timings_$id"
                timings = open((@__DIR__) * "/systems/$name/timings_$id", "r")
            catch e
                @debug "Cannot collect timings for $name / $id"
                push!(cannot_collect, (name, id))
                continue
            end
            lines = readlines(timings)
            if isempty(lines)
                @debug "Cannot collect timings for $name / $id"
                push!(cannot_collect, (name, id))
                continue
            end
            @assert lines[1] == name
            for line in lines[2:end]
                k, v = split(line, ", ")
                runtimes[name][id][Symbol(k)] = parse(Float64, v)
            end
            close(timings)

            #####

            data[name][id] = []
            res = nothing
            try
                @debug "==== Opening /systems/$name/data_$id"
                res = open((@__DIR__) * "/systems/$name/data_$id", "r")
            catch e
                @debug "Cannot collect data for $name / $id"
                push!(cannot_collect, (name, id))
                continue
            end
            lines = readlines(res)
            if isempty(lines)
                @debug "Cannot collect data for $name / $id"
                push!(cannot_collect, (name, id))
                continue
            end
            @assert lines[1] == name
            for line in lines[2:end]
                k, v = split(line, ", ")
                push!(data[name][id], (Symbol(k), v))
            end
            close(res)
        end
    end

    if !isempty(cannot_collect)
        printstyled("(!) Cannot collect timings for:\n", color = :red)
        for (name, id) in cannot_collect
            println("\t$name / $id")
        end
    end

    if args["reparam"]
        ids = map(keywords_to_id, kwargs)
        resulting_md *=
            "|Model|" *
            join(
                vcat("Total time", map(c -> HUMAN_READABLE_CATEGORIES[c], categories)),
                "|",
            ) *
            "|\n"
        resulting_md *=
            "|-----|" * join(["---" for _ in vcat(ids, categories)], "|") * "|\n"
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
            for c in ids
                if !haskey(data[name], c)
                    for _ in 1:length(categories)
                        resulting_md *= " - " * "|"
                    end
                else
                    for (k, v) in data[name][c]
                        if k == :implicit_relations
                            v = parse(Bool, v)
                            v = v ? "yes" : "no"
                            resulting_md *= " $v " * "|"
                        else
                            resulting_md *= " $v " * "|"
                        end
                    end
                end
            end

            resulting_md *= "\n"
        end
    else
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
    timestamp = time_ns()
    args = parse_commandline()
    kwargs = parse_keywords(args["keywords"])
    @debug "Command-line args:"
    for (arg, val) in args
        @debug "$arg  =>  $val"
    end
    @debug "Keywords for `find_identifiable_functions`"
    @debug kwargs

    flag = populate_benchmarks(args, kwargs)

    systems = run_benchmarks(args, kwargs)

    collect_timings(args, kwargs, systems)

    printstyled(
        "Benchmarking finished in $(round((time_ns() - timestamp) / 1e9, digits=2)) s\n",
        color = :light_green,
    )
end

main()
