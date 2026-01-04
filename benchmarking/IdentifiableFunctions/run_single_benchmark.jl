using CpuId, Logging, Pkg, Printf
using Statistics

using StructuralIdentifiability
using StructuralIdentifiability: _runtime_logger, ODE
using StructuralIdentifiability.ParamPunPam

logger = Logging.ConsoleLogger(stdout, Logging.Info)
global_logger(logger)

include("utils.jl")

const data = Dict()

const NUM_RUNS = 5
const FUNCTION_NAME = ARGS[1]
const PROBLEM_NAME = ARGS[2]
const KWARGS = parse_keywords(ARGS[3])[1]
const GLOBAL_ID = keywords_to_global_id(KWARGS)

@info "" FUNCTION_NAME
@info "" PROBLEM_NAME
@info "" KWARGS
@info "" GLOBAL_ID
flush(stdout)
flush(stderr)

# Load the system
path = (@__DIR__) * "/$BENCHMARK_RESULTS/$PROBLEM_NAME/$PROBLEM_NAME.jl"
include(path)

# Compile
const FUNC_SYM = Symbol(FUNCTION_NAME)
ode = @ODEmodel(x1'(t) = a * x1(t) + x2, x2'(t) = b * x2(t) + c * d, y(t) = x1(t))
eval(:($FUNC_SYM($system; ($KWARGS...))))

function process_system()
    @info "Processing $PROBLEM_NAME"
    @info """
    Averaging over $NUM_RUNS runs.
    Using keyword arguments:
    $(typeof(KWARGS))
    $KWARGS
    ID: $GLOBAL_ID"""

    data[PROBLEM_NAME] = Dict{Any, Any}(c => 0.0 for c in ALL_CATEGORIES)
    for _ in 1:NUM_RUNS
        timing = @timed result = eval(:($FUNC_SYM($system; ($KWARGS...))))
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
    return
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
    return open((@__DIR__) * "/$BENCHMARK_RESULTS/$PROBLEM_NAME/$filename", "w") do io
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
    open((@__DIR__) * "/$BENCHMARK_RESULTS/$PROBLEM_NAME/$filename", "w") do io
        write(io, "$PROBLEM_NAME\n")
    end
    for cat in ID_DATA_CATEGORIES
        if !haskey(data[PROBLEM_NAME], cat)
            continue
        end
        if cat === :something_important
            # make a separate file for it
            filename_cat = generic_filename(cat, GLOBAL_ID)
            open((@__DIR__) * "/$BENCHMARK_RESULTS/$PROBLEM_NAME/$filename_cat", "w") do io
                # print something
            end
            continue
        end
        # otherwise, print in the data file
        open((@__DIR__) * "/$BENCHMARK_RESULTS/$PROBLEM_NAME/$filename", "a+") do io
            write(io, "$cat, ")
            write(io, string(data[PROBLEM_NAME][cat]))
            write(io, "\n")
        end
    end
    return
end

process_system()
dump_timings()
dump_results()
