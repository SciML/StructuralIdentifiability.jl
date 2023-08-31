using CpuId, Logging, Pkg, Printf
using Statistics

using StructuralIdentifiability
using StructuralIdentifiability: _runtime_logger, ODE
using StructuralIdentifiability.ParamPunPam

logger = Logging.SimpleLogger(stdout, Logging.Warn)
global_logger(logger)

include("utils.jl")

runtimes = Dict()
data = Dict()
results = Dict()

NUM_RUNS = 3
NAME = ARGS[1]
KWARGS = parse_keywords(ARGS[2])[1]
ID = keywords_to_id(KWARGS)

# Load the system
path = (@__DIR__) * "/systems/$NAME/$NAME.jl"
include(path)

# Compile
ode = @ODEmodel(
    x1'(t) = (-b * c * x1(t) - b * x1(t) * x4(t) + 1) // (c + x4(t)),
    x2'(t) = alpha * x1(t) - beta * x2(t),
    x3'(t) = gama * x2(t) - delta * x3(t),
    x4'(t) = (gama * sigma * x2(t) * x4(t) - delta * sigma * x3(t) * x4(t)) // x3(t),
    y(t) = x1(t)
)
find_identifiable_functions(ode; KWARGS...)

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

function process_system()
    @info "Processing $NAME"
    @info """
    Averaging over $NUM_RUNS runs.
    Using keyword arguments:
    $(typeof(KWARGS))
    $KWARGS
    ID: $ID"""

    runtimes[NAME] = Dict(c => 0.0 for c in ALL_CATEGORIES)
    data[NAME] = Dict{Any, Any}(c => [] for c in ID_DATA_CATEGORIES)
    for _ in 1:NUM_RUNS
        funcs = find_identifiable_functions(system; KWARGS...)
        results[NAME] = funcs
        @info "Identifiable functions are" funcs
        for cat in ID_TIME_CATEGORIES
            runtimes[NAME][cat] += StructuralIdentifiability._runtime_logger[cat]
        end
        for cat in ID_DATA_CATEGORIES
            data[NAME][cat] = deepcopy(StructuralIdentifiability._runtime_logger[cat])
        end
    end
    for k in keys(runtimes[NAME])
        runtimes[NAME][k] = runtimes[NAME][k] / NUM_RUNS
    end
end

function dump_timings()
    timings = ""
    timings *= "$NAME\n"

    for (name, times) in runtimes
        for c in ALL_CATEGORIES
            timings *= "$c, "
            timings *= string(times[c]) * "\n"
        end
    end

    open((@__DIR__) * "/systems/$NAME/timings_$ID", "w") do io
        write(io, timings)
    end
end

function dump_results()
    # dump identifiabile functions
    open((@__DIR__) * "/systems/$NAME/id_funcs_$ID", "w") do io
        if haskey(results, NAME)
            funcs_str = "[" * join(map(repr, results[NAME]), ",\n") * "]"
            println(io, funcs_str)
        end
    end
    # dump everything else
    for cat in ID_DATA_CATEGORIES
        open((@__DIR__) * "/systems/$NAME/$(cat)", "w") do io
            if haskey(data, NAME)
                if cat === :id_certain_factors
                    factors = data[NAME][cat]
                    number = length(factors)
                    med = median(map(length, factors))
                    R = parent(factors[1][1])
                    factors_str = map(s -> map(repr, s), factors)
                    factors_str = map(f -> join(f, ", "), factors_str)
                    factors_str = join(factors_str, "\n")
                    println(io, "Factors in $R")
                    println(io, "Factored $number polynomials in total")
                    println(io, "The median number of factors is $med")
                    println(io, "\n============================\n")
                    println(io, factors_str)
                else
                    @warn "Skipping printing data: $cat"
                end
            end
        end
    end
end

process_system()
dump_timings()
dump_results()
