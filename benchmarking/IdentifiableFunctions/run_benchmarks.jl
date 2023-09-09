#=
This script runs benchmarks and collects useful statistics.

To benchmark a particular function `StructuralIdentifiability.XYZ`, execute the
following command in your favorite terminal from this directory:

```
$ julia run_benchmarks.jl XYZ
```

`XYZ` can be one of the following:
- `find_identifiable_functions`,
- `reparametrize_global`.

Going into details, this will:
1. Load all benchmark models from benchmarks.jl and create a directory
   `benchmark_results`.
2. Run the function `StructuralIdentifiability.XYZ`, collect useful runtime
   statistics, and write them to `benchmark_results`.
3. Read the statistics from `benchmark_results` and produce a table with
   results.

Below in the file there is a description of some potentially useful command line
arguments.
A more involved way of running benchmarks could be the following:

```
julia run_benchmarks.jl XYZ --target="id_total, id_are_polynomial" --timeout=3600 --keywords="(strategy=(:hybrid,),with_states=true); (strategy=(:normalforms,3),with_states=true)"
```
=#

using ArgParse
using CpuId, Logging, Pkg, Printf
using Base.Threads
using Distributed
using Dates
using ProgressMeter

using StructuralIdentifiability
using StructuralIdentifiability: _runtime_logger, ODE
using StructuralIdentifiability.ParamPunPam

global_logger(Logging.ConsoleLogger(stdout, Logging.Warn))
include("benchmarks.jl")
global_logger(Logging.ConsoleLogger(stdout, Logging.Info))

const _progressbar_color = :light_green
const _progressbar_value_color = :light_green
progressbar_enabled() =
    Logging.Info <= Logging.min_enabled_level(current_logger()) < Logging.Warn

include("utils.jl")

const BENCHMARK_TABLE = "benchmark_result.md"

function parse_commandline()
    s = ArgParseSettings()
    #! format: off
    @add_arg_table s begin
        "function"
            help = "The function to benchmark."
            arg_type = String
            required = true
        "--keywords"
            help = """
            Keyword arguments to be passed to `function`. 
            Semicolon-separated list of named tuples."""
            arg_type = String
            default = ""
            # default = "(strategy=(:gb, ),); (strategy=(:gb, ),with_states=true); (strategy=(:normalforms, 2),); (strategy=(:normalforms, 2),with_states=true); (strategy=(:normalforms, 3),); (strategy=(:normalforms, 3),with_states=true); (strategy=(:hybrid, ),); (strategy=(:hybrid, ),with_states=true)"
        "--target"
            help = """
            Statistics to be displayed in the table.
            Comma-separated list of entities.
            NOTE: each of these must be present in the list of
            `ALL_POSSIBLE_CATEGORIES` in `utils.jl`"""
            arg_type = String
            default = "id_total"
        "--timeout"
            help = "Timeout, s."
            arg_type = Int
            default = 300
        "--workers"
            help = "The number of available worker processes."
            arg_type = Int
            default = 4
        "--skip"
            help = "Skip specified benchmark models."
            arg_type = Vector{String}
            default = ["NFkB"]
        "--augment"
            help = "Augment the benchmark dataset with similar models."
            arg_type = Bool
            default = false
        "--regen"
            help = "Re-generate the folder with benchmarks from scratch."
            arg_type = Bool
            default = false
    end
    #! format: on

    parse_args(s)
end

function populate_benchmarks(args, kwargs)
    regen = args["regen"]
    dir_present = isdir((@__DIR__) * "/$BENCHMARK_RESULTS/")
    dir_present && !regen && return false
    @debug "Re-generating the benchmarks folder"
    try
        if isdir((@__DIR__) * "/$BENCHMARK_RESULTS/")
            rm((@__DIR__) * "/$BENCHMARK_RESULTS/", recursive = true, force = true)
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
        mkpath((@__DIR__) * "/$BENCHMARK_RESULTS/$name/")
        fd = open((@__DIR__) * "/$BENCHMARK_RESULTS/$name/$name.jl", "w")
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

function run_benchmarks(args, kwargs)
    to_skip = args["skip"]
    timeout = args["timeout"]
    @assert timeout > 0
    function_name = args["function"]
    @assert function_name in ("find_identifiable_functions", "reparametrize_global")
    nworkers = args["workers"]
    @assert nworkers > 0

    dirnames = first(walkdir((@__DIR__) * "/$BENCHMARK_RESULTS/"))[2]
    to_run_names = setdiff(dirnames, to_skip)[1:3]
    to_run_indices = collect(1:length(to_run_names))

    @info """
    Benchmarking `$function_name`."""
    @info """
    Passing these keyword arguments to `$function_name`:
    \t$(join(map(string, kwargs), "\n\t"))"""
    @info """
    Number of benchmark systems: $(length(to_run_indices))
    Workers: $(nworkers)
    Timeout: $timeout seconds"""
    @info """
    Benchmark systems:
    $to_run_names"""

    seconds_passed(from_t) = round((time_ns() - from_t) / 1e9, digits = 2)

    queue = [
        (problem_id = problem, function_kwargs = kw) for kw in kwargs for
        problem in to_run_indices
    ]
    processes = []
    running = []
    errored = []

    generate_showvalues(processes) =
        () -> [(
            :Active,
            join(
                map(
                    proc ->
                        string(proc.problem_name) * " / " * string(proc.global_run_id),
                    filter(proc -> process_running(proc.julia_process), processes),
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
        if !isempty(queue) && length(running) < nworkers
            next!(
                prog,
                showvalues = generate_showvalues(processes),
                step = 0,
                valuecolor = _progressbar_value_color,
                # spinner = "⌜⌝⌟⌞",
            )
            task = pop!(queue)
            function_kwargs = task.function_kwargs
            problem_id = task.problem_id
            problem_name = to_run_names[problem_id]
            global_id = keywords_to_global_id(function_kwargs)
            @debug "Running $problem_name / $global_id"
            logfn = global_id === Symbol("") ? "logs" : "logs_$global_id"
            errfn = global_id === Symbol("") ? "errors" : "errors_$global_id"
            logs = open((@__DIR__) * "/$BENCHMARK_RESULTS/$problem_name/$logfn", "w")
            errs = open((@__DIR__) * "/$BENCHMARK_RESULTS/$problem_name/$errfn", "w")
            cmd = Cmd([
                "julia",
                (@__DIR__) * "/run_single_benchmark.jl",
                "$function_name",
                "$problem_name",
                "$function_kwargs",
            ])
            cmd = Cmd(cmd, ignorestatus = true, detach = false, env = copy(ENV))
            proc = run(pipeline(cmd, stdout = logs, stderr = errs), wait = false)
            push!(
                processes,
                (
                    problem_id = problem_id,
                    problem_name = problem_name,
                    function_name = function_name,
                    function_kwargs = function_kwargs,
                    julia_process = proc,
                    start_time = time_ns(),
                    global_run_id = global_id,
                    logfile = logs,
                    errfile = errs,
                ),
            )
            push!(running, processes[end])
        end

        sleep(0.2)
        to_be_removed = []
        for i in 1:length(running)
            proc = running[i]
            if process_exited(proc.julia_process)
                push!(to_be_removed, i)
                if proc.julia_process.exitcode != 0
                    push!(errored, proc)
                end
                close(proc.logfile)
                close(proc.errfile)
                kw = proc.function_kwargs
                start_time = proc.start_time
                next!(
                    prog,
                    showvalues = generate_showvalues(processes),
                    valuecolor = _progressbar_value_color,
                )
                @debug "Yielded $(proc.problem_name) / $(kw) after $(seconds_passed(start_time)) seconds"
            end
            if process_running(proc.julia_process)
                start_time = proc.start_time
                if seconds_passed(start_time) > timeout
                    push!(to_be_removed, i)
                    kill(proc.julia_process)
                    close(proc.logfile)
                    close(proc.errfile)
                    kw = proc.function_kwargs
                    next!(
                        prog,
                        showvalues = generate_showvalues(processes),
                        valuecolor = _progressbar_value_color,
                    )
                    @debug "Timed-out $(proc.problem_name) / $(kw) after $(seconds_passed(start_time)) seconds"
                end
            end
        end
        deleteat!(running, to_be_removed)
        if isempty(queue)
            @debug "All benchmarks finished"
            break
        end
    end
    finish!(prog)

    if !isempty(errored)
        printstyled("(!) Maybe errored:\n", color = :red)
        for proc in errored
            println("\t$(proc.problem_name) / $(proc.global_run_id)")
        end
    end

    to_run_names
end

function collect_timings(args, kwargs, names; content = :compare)
    function_name = args["function"]
    targets = map(Symbol, map(strip, split(args["target"], ",")))
    @assert all(target -> target in ALL_CATEGORIES, targets)
    @assert length(targets) > 0
    @info """
    Collecting benchmark results for `$function_name`.

    Keyword arguments of interest:
    \t$(join(map(string, kwargs), "\n\t"))

    Statistics of interest:
    \t$(join(map(string, targets), "\n\t"))
    """

    cannot_collect = []
    names = sort(names)
    kwids = map(keywords_to_global_id, kwargs)

    # Collect timings and data from directory BENCHMARK_RESULTS.
    data = Dict()
    for name in names
        @debug "==== Reading $name"
        data[name] = Dict()
        for kwid in kwids
            timings_file = nothing
            #####
            data[name][kwid] = Dict()
            timingsfn = timings_filename(kwid)
            try
                @debug "==== Opening /$BENCHMARK_RESULTS/$name/$timingsfn"
                timings_file =
                    open((@__DIR__) * "/$BENCHMARK_RESULTS/$name/$timingsfn", "r")
            catch e
                @debug "Cannot collect timings for $name / $kwid"
                push!(cannot_collect, (name, kwid))
                continue
            end
            lines = readlines(timings_file)
            if isempty(lines)
                @debug "Cannot collect timings for $name / $kwid"
                push!(cannot_collect, (name, kwid))
                continue
            end
            @assert lines[1] == name
            for line in lines[2:end]
                k, v = split(line, ", ")
                data[name][kwid][Symbol(k)] = parse(Float64, v)
            end
            close(timings_file)
            #####
            datafn = data_filename(kwid)
            data_file = nothing
            try
                @debug "==== Opening /$BENCHMARK_RESULTS/$name/$datafn"
                data_file = open((@__DIR__) * "/$BENCHMARK_RESULTS/$name/$datafn", "r")
            catch e
                @debug "Cannot collect data for $name / $kwid"
                push!(cannot_collect, (name, kwid))
                continue
            end
            lines = readlines(data_file)
            if isempty(lines)
                @debug "Cannot collect data for $name / $kwid"
                push!(cannot_collect, (name, kwid))
                continue
            end
            @assert lines[1] == name
            for line in lines[2:end]
                k, v = map(strip, split(line, ","))
                data[name][kwid][Symbol(k)] = v
            end
            close(data_file)
        end
    end

    if !isempty(cannot_collect)
        printstyled("(!) Cannot collect timings for:\n", color = :red)
        for (name, kwid) in cannot_collect
            println("\t$name / $kwid")
        end
    end

    # Print the table to BENCHMARK_TABLE.
    resulting_md = ""
    resulting_md *= """
    ## Benchmark results

    $(now())

    Function: `$(args["function"])`

    Workers: $(args["workers"])

    Timeout: $(args["timeout"]) s

    **All timings in seconds.**

    """

    makecolname(kw, target) =
        kw === Symbol("") ? HUMAN_READABLE_CATEGORIES[target] :
        Symbol(Symbol(kw), Symbol("/"), HUMAN_READABLE_CATEGORIES[target])
    columns = [makecolname(kwid, target) for kwid in kwids for target in targets]
    resulting_md *= "|Model|" * join(map(string, columns), "|") * "|\n"
    resulting_md *= "|-----|" * join(["---" for _ in columns], "|") * "|\n"
    for name in names
        model_data = data[name]
        resulting_md *= "|$name|"
        for kwid in kwids
            if !haskey(model_data, kwid)
                resulting_md *= (" - " * "|")^length(columns)
                continue
            end
            for target in targets
                if !haskey(model_data[kwid], target)
                    resulting_md *= " - " * "|"
                else
                    formatting_style = CATEGORY_FORMAT[target]
                    resulting_md *= formatting_style(model_data[kwid][target]) * "|"
                end
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
        if deps[uid].version !== nothing
            resulting_md *= "* $s : $(deps[uid].version)\n"
        end
    end

    open((@__DIR__) * "/$BENCHMARK_TABLE", "w") do io
        write(io, resulting_md)
    end
end

function main()
    timestamp = time_ns()
    args = parse_commandline()
    @debug "Command-line args:"
    for (arg, val) in args
        @debug "$arg  =>  $val"
    end
    kwargs = parse_keywords(args["keywords"])
    @debug "Keywords for `$(args["function"])`"
    @debug kwargs
    flag = populate_benchmarks(args, kwargs)
    problems = run_benchmarks(args, kwargs)
    printstyled(
        """
        Benchmarking had finished in $(round((time_ns() - timestamp) / 1e9, digits=2)) seconds.
        Results are written to /$BENCHMARK_RESULTS
        """,
        color = :light_green,
    )
    collect_timings(args, kwargs, problems)
    printstyled(
        """
        Table with results is written to /$BENCHMARK_TABLE
        """,
        color = :light_green,
    )
end

main()
