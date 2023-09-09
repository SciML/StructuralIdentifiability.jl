using CpuId, Logging, Pkg, Printf
using Statistics

using StructuralIdentifiability
using StructuralIdentifiability: _runtime_logger, ODE
using StructuralIdentifiability.ParamPunPam

logger = Logging.ConsoleLogger(stdout, Logging.Info)
global_logger(logger)

include("utils.jl")

const data = Dict()

const NUM_RUNS = 1
const FUNCTION_NAME = ARGS[1]
const PROBLEM_NAME = ARGS[2]
const KWARGS = parse_keywords(ARGS[3])[1]
const GLOBAL_ID = keywords_to_global_id(KWARGS)

# Load the system
path = (@__DIR__) * "/$BENCHMARK_RESULTS/$PROBLEM_NAME/$PROBLEM_NAME.jl"
include(path)

# Compile
# ode = @ODEmodel(
#     x1'(t) = (-b * c * x1(t) - b * x1(t) * x4(t) + 1) // (c + x4(t)),
#     x2'(t) = alpha * x1(t) - beta * x2(t),
#     x3'(t) = gama * x2(t) - delta * x3(t),
#     x4'(t) = (gama * sigma * x2(t) * x4(t) - delta * sigma * x3(t) * x4(t)) // x3(t),
#     y(t) = x1(t)
# )
# find_identifiable_functions(ode; KWARGS...)

macro invoke_function(func, args, kwargs)
    func_sym = Symbol(func)
    esc(:($func_sym($args, $(kwargs)...)))
end

function process_system()
    @info "Processing $PROBLEM_NAME"
    @info """
    Averaging over $NUM_RUNS runs.
    Using keyword arguments:
    $(typeof(KWARGS))
    $KWARGS
    ID: $GLOBAL_ID"""

    func_sym = Symbol(FUNCTION_NAME)
    data[PROBLEM_NAME] = Dict{Any, Any}(c => 0.0 for c in ALL_CATEGORIES)
    for _ in 1:NUM_RUNS
        timing = @timed result = eval(:($func_sym($system, $(KWARGS...))))
        data[PROBLEM_NAME][:return_value] = result
        @info "Result is" result
        for cat in ID_TIME_CATEGORIES
            if haskey(StructuralIdentifiability._runtime_logger, cat)
                data[PROBLEM_NAME][cat] = StructuralIdentifiability._runtime_logger[cat]
            end
        end
        for cat in ID_DATA_CATEGORIES
            if haskey(StructuralIdentifiability._runtime_logger, cat)
                data[PROBLEM_NAME][cat] =
                    deepcopy(StructuralIdentifiability._runtime_logger[cat])
            end
        end
        data[PROBLEM_NAME][:id_total] = timing.time
    end
    for cat in ID_TIME_CATEGORIES
        if haskey(data[PROBLEM_NAME], cat)
            data[PROBLEM_NAME][cat] = data[PROBLEM_NAME][cat] / NUM_RUNS
        end
    end
end

function dump_timings()
    timings = ""
    timings *= "$PROBLEM_NAME\n"
    for (key, model_data) in data
        for c in ID_TIME_CATEGORIES
            timings *= "$c, "
            timings *= string(model_data[c]) * "\n"
        end
    end
    filename = timings_filename(GLOBAL_ID)
    open((@__DIR__) * "/$BENCHMARK_RESULTS/$PROBLEM_NAME/$filename", "w") do io
        write(io, timings)
    end
end

function dump_results()
    filename = result_filename(GLOBAL_ID)
    open((@__DIR__) * "/$BENCHMARK_RESULTS/$PROBLEM_NAME/$filename", "w") do io
        if haskey(data, PROBLEM_NAME)
            println(io, data[PROBLEM_NAME][:return_value])
        end
    end
    filename = data_filename(GLOBAL_ID)
    for cat in ID_DATA_CATEGORIES
        if cat === :something_important
            # make a separate file for it
            filename_cat = generic_filename(cat, GLOBAL_ID)
            open((@__DIR__) * "/$BENCHMARK_RESULTS/$NAME/$filename_cat", "w") do io
                # print something
            end
            continue
        end
        # otherwise, print in the data file
        open((@__DIR__) * "/$BENCHMARK_RESULTS/$NAME/$filename", "w+") do io
            write(io, string(model_data[cat]))
        end
    end
end

process_system()
dump_timings()
dump_results()
