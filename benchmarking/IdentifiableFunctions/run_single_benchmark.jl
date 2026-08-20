using Logging
using StructuralIdentifiability

global_logger(Logging.ConsoleLogger(stdout, Logging.Info))

include("utils.jl")

const NUM_RUNS = 1
const FUNCTION_NAME = ARGS[1]
const PROBLEM_NAME = ARGS[2]
const KWARGS = only(parse_keywords(ARGS[3]))
const GLOBAL_ID = keywords_to_global_id(KWARGS)
const RESULT_DIR = joinpath(@__DIR__, BENCHMARK_RESULTS, PROBLEM_NAME)

@info "" FUNCTION_NAME
@info "" PROBLEM_NAME
@info "" KWARGS
@info "" GLOBAL_ID
flush(stdout)
flush(stderr)

include(joinpath(RESULT_DIR, "$PROBLEM_NAME.jl"))

const RESOLVED_KWARGS = if haskey(KWARGS, :cmp) && KWARGS.cmp isa Symbol
    merge(KWARGS, (cmp = benchmark_comparator(KWARGS.cmp, system),))
else
    KWARGS
end
const BENCHMARK_FUNCTION = getproperty(StructuralIdentifiability, Symbol(FUNCTION_NAME))

# Compile before collecting timings.
BENCHMARK_FUNCTION(system; RESOLVED_KWARGS...)

function process_system()
    @info "Processing $PROBLEM_NAME"
    @info """
    Averaging over $NUM_RUNS runs.
    Using keyword arguments:
    $RESOLVED_KWARGS
    ID: $GLOBAL_ID"""

    model_data = Dict{Symbol, Any}(category => 0.0 for category in ID_TIME_CATEGORIES)
    for _ in 1:NUM_RUNS
        timing = @timed BENCHMARK_FUNCTION(system; RESOLVED_KWARGS...)
        model_data[:return_value] = timing.value
        @info "Result is" timing.value

        for category in ID_TIME_CATEGORIES
            category === :id_total && continue
            if haskey(StructuralIdentifiability._runtime_logger, category)
                model_data[category] += StructuralIdentifiability._runtime_logger[category]
            end
        end
        for category in ID_DATA_CATEGORIES
            if haskey(StructuralIdentifiability._runtime_logger, category)
                model_data[category] =
                    deepcopy(StructuralIdentifiability._runtime_logger[category])
            end
        end
        model_data[:id_total] += timing.time
    end

    for category in ID_TIME_CATEGORIES
        model_data[category] /= NUM_RUNS
    end
    return model_data
end

function dump_timings(model_data)
    filename = joinpath(RESULT_DIR, timings_filename(GLOBAL_ID))
    return open(filename, "w") do io
        println(io, PROBLEM_NAME)
        for category in ID_TIME_CATEGORIES
            println(io, "$category, $(model_data[category])")
        end
    end
end

function dump_results(model_data)
    result_path = joinpath(RESULT_DIR, result_filename(GLOBAL_ID))
    open(result_path, "w") do io
        println(io, model_data[:return_value])
    end

    data_path = joinpath(RESULT_DIR, data_filename(GLOBAL_ID))
    return open(data_path, "w") do io
        println(io, PROBLEM_NAME)
        for category in ID_DATA_CATEGORIES
            haskey(model_data, category) || continue
            println(io, "$category, $(model_data[category])")
        end
    end
end

const MODEL_DATA = process_system()
dump_timings(MODEL_DATA)
dump_results(MODEL_DATA)
