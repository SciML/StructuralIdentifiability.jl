#=
Run identifiable-function benchmarks and compare the returned generators.

From this directory, for example:

```sh
julia run_benchmarks.jl find_identifiable_functions \
    --timeout=3600 \
    --memorylimit=8192 \
    --keywords="(with_states=true, cmp=:default); (with_states=true, cmp=:cmp_lie)"
```

Generated model files and results are written below `results/`; subprocess logs
are written below `logs/`. Both directories are ignored by Git. The summarized
comparison is written to `benchmark_result.md`.
=#

using ArgParse
using CpuId, Dates, Logging, Pkg, ProgressMeter
using StructuralIdentifiability

global_logger(Logging.ConsoleLogger(stdout, Logging.Warn))
include("../benchmarks.jl")
global_logger(Logging.ConsoleLogger(stdout, Logging.Info))

const _progressbar_color = :light_green
const _progressbar_value_color = :light_green
progressbar_enabled() =
    Logging.Info <= Logging.min_enabled_level(current_logger()) < Logging.Warn

include("utils.jl")

const RESULTS_DIR = joinpath(@__DIR__, BENCHMARK_RESULTS)
const LOGS_DIR = joinpath(@__DIR__, BENCHMARK_LOGS)
const BENCHMARK_TABLE = joinpath(@__DIR__, "benchmark_result.md")

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
        "--timeout"
            help = "Timeout, s."
            arg_type = Int
            default = 3600
        "--memorylimit"
            help = "Optional memory limit in MiB for each benchmark subprocess. Use 0 to disable."
            arg_type = Int
            default = 0
        "--workers"
            help = "The number of available worker processes."
            arg_type = Int
            default = 4
        "--skip"
            help = "Skip specified benchmark models."
            arg_type = Vector{Symbol}
            default = [
                :NFkB,
                :Covid2,
                :Lincomp1,
                :Lincomp2,
                Symbol("cLV1 (1o)"),
                Symbol("cLV1 (2o)"),
                :TumorHu2019,
                :TumorPillis2007,
            ]
        "--models"
            help = """
            Run specified benchmark models. 
            A comma-separated list of benchmark IDs.
            Leave empty for selecting all models."""
            arg_type = String
            default = ""
        "--regen"
            help = "Re-generate the folder with benchmarks from scratch."
            arg_type = Bool
            default = false
        "--tableonly"
            help = """
            Do not run benchmarks. 
            Just construct the table from the existing directory."""
            arg_type = Bool
            default = false
    end
    #! format: on

    return parse_args(s)
end

function selected_benchmark_names(args)
    requested_ids = filter(!isempty, strip.(split(args["models"], ",")))
    names = if isempty(requested_ids)
        filter(name -> isdir(joinpath(RESULTS_DIR, name)), readdir(RESULTS_DIR))
    else
        [benchmarks[Symbol(id)][:name] for id in requested_ids]
    end
    skipped_names = [
        benchmarks[id][:name] for id in args["skip"] if haskey(benchmarks, id)
    ]
    return sort(setdiff(names, skipped_names))
end

function populate_benchmarks(args)
    results_present = isdir(RESULTS_DIR)
    results_present && !args["regen"] && return
    results_present && rm(RESULTS_DIR; recursive = true, force = true)
    mkpath(RESULTS_DIR)

    @debug "Re-generating the benchmarks folder"
    progress = Progress(
        length(benchmarks),
        "Generating benchmarks";
        dt = 0.1,
        enabled = progressbar_enabled(),
        color = _progressbar_color,
    )
    for (_, benchmark) in benchmarks
        next!(progress)
        name = benchmark[:name]
        system = benchmark[:ode]
        @debug "Generating $name"
        model_dir = joinpath(RESULTS_DIR, name)
        mkpath(model_dir)
        open(joinpath(model_dir, "$name.jl"), "w") do io
            println(io, "# $name")
            println(io, "#! format: off")
            println(io, "using StructuralIdentifiability\n")
            equations = join(
                map(line -> "\t" * line, split(repr(system), "\n")[1:(end - 1)]),
                ",\n",
            )
            println(io, "system = @ODEmodel(\n$equations\n)")
        end
    end
    finish!(progress)
    return
end

function benchmark_command(function_name, problem_name, function_kwargs, memory_limit)
    project_dir = dirname(Base.active_project())
    script = joinpath(@__DIR__, "run_single_benchmark.jl")
    command = `$(Base.julia_cmd()) --project=$project_dir $script $function_name $problem_name $(string(function_kwargs))`
    command = Cmd(command; ignorestatus = true, detach = false, env = copy(ENV))
    memory_limit == 0 && return command

    Sys.isunix() || throw(ArgumentError("--memorylimit is supported only on Unix"))
    shell_command = string(
        "ulimit -v ",
        memory_limit * 1024,
        "; exec ",
        join(Base.shell_escape.(command.exec), " "),
    )
    shell = Cmd(["bash", "-lc", shell_command])
    return Cmd(shell; ignorestatus = true, detach = false, env = copy(ENV))
end

function run_benchmarks(args, kwargs)
    timeout = args["timeout"]
    @assert timeout > 0
    memorylimit = args["memorylimit"]
    @assert memorylimit >= 0
    function_name = args["function"]
    @assert function_name in ("find_identifiable_functions", "reparametrize_global")
    nworkers = args["workers"]
    @assert nworkers > 0

    to_run_names = selected_benchmark_names(args)
    to_run_indices = eachindex(to_run_names)

    @info """
    Benchmarking `$function_name`."""
    @info """
    Passing these keyword arguments to `$function_name`:
    $(join(map(string, kwargs), "\n\t"))"""
    @info """
    Number of benchmark systems: $(length(to_run_indices))
    Workers: $(nworkers)
    Timeout: $timeout seconds"""
    @info """
    Benchmark systems:
    $to_run_names"""

    seconds_passed(from_t) = round((time_ns() - from_t) / 1.0e9, digits = 2)

    queue = [
        (problem_id = problem, function_kwargs = kw) for kw in kwargs for
            problem in to_run_indices
    ]
    processes = []
    running = []
    errored = []

    generate_showvalues(processes) =
        () -> [
        (
            :Active,
            join(
                map(
                    proc ->
                    string(proc.problem_name) * " / " * string(proc.global_run_id),
                    filter(proc -> process_running(proc.julia_process), processes),
                ),
                ", ",
            ),
        ),
    ]

    prog = Progress(
        length(queue),
        "Running benchmarks",
        dt = 0.3,
        enabled = progressbar_enabled(),
        color = _progressbar_color,
    )
    while true
        if !isempty(queue) && length(running) < nworkers
            task = pop!(queue)
            function_kwargs = task.function_kwargs
            problem_id = task.problem_id
            problem_name = to_run_names[problem_id]
            global_id = keywords_to_global_id(function_kwargs)
            @debug "Running $problem_name / $global_id. Kwargs:\n$function_kwargs"
            log_dir = joinpath(LOGS_DIR, problem_name)
            mkpath(log_dir)
            log_path = joinpath(log_dir, generic_filename("run", global_id) * ".log")
            logs = open(log_path, "w")
            command = benchmark_command(
                function_name,
                problem_name,
                function_kwargs,
                memorylimit,
            )
            proc = run(pipeline(command; stdout = logs, stderr = logs); wait = false)
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
                ),
            )
            push!(running, processes[end])
            next!(
                prog,
                showvalues = generate_showvalues(processes),
                step = 0,
                valuecolor = _progressbar_value_color,
            )
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
                    wait(proc.julia_process)
                    close(proc.logfile)
                    push!(errored, proc)
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
        if isempty(queue) && isempty(running)
            @debug "All benchmarks finished"
            break
        end
    end
    finish!(prog)

    if !isempty(errored)
        printstyled("Failed or timed out:\n", color = :red)
        for proc in errored
            println("\t$(proc.problem_name) / $(proc.global_run_id)")
        end
    end

    return to_run_names
end

function compact_result_string(content)
    content = strip(content)
    isempty(content) && return "-"
    first_bracket = findfirst(==('['), content)
    if first_bracket !== nothing
        content = content[first_bracket:end]
    end
    content = replace(content, "\n" => " ")
    content = replace(content, r"\s+" => " ")
    content = replace(content, "|" => "\\|")
    return strip(content)
end

function display_keyword_label(keyword)
    isempty(keyword) && return "default"
    if haskey(keyword, :cmp)
        cmp_label = string(keyword.cmp)
        if haskey(keyword, :with_states) && keyword.with_states
            return "funcs ($cmp_label)"
        end
        return string(cmp_label)
    end
    return string(keywords_to_global_id(keyword))
end

function collect_results(args, kwargs, names)
    function_name = args["function"]
    @info """
    Collecting benchmark results for `$function_name`.

    Keyword arguments of interest:
    $(join(map(string, kwargs), "\n\t"))
    """

    cannot_collect = []
    names = sort(names)
    kwids = map(keywords_to_global_id, kwargs)
    labels = map(display_keyword_label, kwargs)
    memory_limit = args["memorylimit"]
    memory_description = memory_limit == 0 ? "disabled" : "$memory_limit MiB"

    data = Dict{String, Dict{Symbol, String}}()
    for name in names
        model_data = Dict{Symbol, String}()
        data[name] = model_data
        for id in kwids
            result_path = joinpath(RESULTS_DIR, name, result_filename(id))
            if isfile(result_path)
                model_data[id] = compact_result_string(read(result_path, String))
            else
                @debug "Cannot collect result for $name / $id"
                push!(cannot_collect, (name, id))
                model_data[id] = "-"
            end
        end
    end

    if !isempty(cannot_collect)
        printstyled("(!) Cannot collect results for:\n", color = :red)
        for (name, kwid) in cannot_collect
            println("\t$name / $kwid")
        end
    end

    resulting_md = ""
    resulting_md *= """
    ## Benchmark results

    $(now())

    - Benchmarked function: `$(args["function"])`
    - Workers: $(args["workers"])
    - Timeout: $(args["timeout"]) s
    - Memory limit: $memory_description

    """

    columns = ["Model"; labels]
    resulting_md *= "|" * join(columns, "|") * "|\n"
    resulting_md *= "|" * join(["---" for _ in columns], "|") * "|\n"
    for name in names
        resulting_md *= "|$name|"
        for kwid in kwids
            resulting_md *= data[name][kwid] * "|"
        end
        resulting_md *= "\n"
    end

    resulting_md *= "\n*Benchmarking environment:*\n\n"
    resulting_md *= "* Total RAM (GiB): $(div(Sys.total_memory(), 2^30))\n"
    resulting_md *= "* Processor: $(cpubrand())\n"
    resulting_md *= "* Julia version: $(VERSION)\n\n"
    resulting_md *= "Versions of the dependencies:\n\n"

    deps = Pkg.dependencies()
    dependencies = sort(collect(Pkg.project().dependencies); by = first)
    for (name, uuid) in dependencies
        version = deps[uuid].version
        version === nothing || (resulting_md *= "* $name : $version\n")
    end

    return open(BENCHMARK_TABLE, "w") do io
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
    if args["tableonly"]
        problems = selected_benchmark_names(args)
    else
        populate_benchmarks(args)
        problems = run_benchmarks(args, kwargs)
        printstyled(
            """
            Benchmarking had finished in $(round((time_ns() - timestamp) / 1.0e9, digits = 2)) seconds.
            Results are written to $RESULTS_DIR
            """,
            color = :light_green,
        )
    end
    collect_results(args, kwargs, problems)
    return printstyled(
        """
        Table with results is written to $BENCHMARK_TABLE
        """,
        color = :light_green,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
