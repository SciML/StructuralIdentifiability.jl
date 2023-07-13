using CpuId, Logging, Pkg, Printf

using StructuralIdentifiability
using StructuralIdentifiability: _runtime_logger, ODE
using StructuralIdentifiability.ParamPunPam

logger = Logging.SimpleLogger(stdout, Logging.Warn)
global_logger(logger)

runtimes = Dict()
results = Dict()

ID_TIME_CATEGORIES = [
    :id_io_time,
    :id_global_time,
    :id_ideal_time,
    :id_filter_time,
    :id_inclusion_check,
    :id_groebner_time,
    :id_total,
]
GB_TIME_CATEGORIES = [
    :gb_discover_shape,
    :gb_discover_degrees,
    :gb_interpolate,
    :gb_recover_coeffs,
    :gb_npoints,
]
ALL_CATEGORIES = union(ID_TIME_CATEGORIES, GB_TIME_CATEGORIES)
NUM_RUNS = 1
NAME = ARGS[1]

# Load the system
path = (@__DIR__) * "/systems/$NAME/$NAME.jl"
include(path)

# Compile
ode = @ODEmodel(
    x1'(t) = (-b*c*x1(t) - b*x1(t)*x4(t) + 1)//(c + x4(t)),
    x2'(t) = alpha*x1(t) - beta*x2(t),
    x3'(t) = gama*x2(t) - delta*x3(t),
    x4'(t) = (gama*sigma*x2(t)*x4(t) - delta*sigma*x3(t)*x4(t))//x3(t),
    y(t) = x1(t)
)
find_identifiable_functions(ode)
for cat in ID_TIME_CATEGORIES
    StructuralIdentifiability._runtime_logger[cat] = 0.0
end

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

function process_system()
    @info "Processing $NAME"
    runtimes[NAME] = Dict(c => 0.0 for c in ALL_CATEGORIES)
    for _ in 1:NUM_RUNS
        funcs = find_identifiable_functions(system)
        results[NAME] = funcs
        @info "Identifiable functions are" funcs
        for cat in ID_TIME_CATEGORIES
            runtimes[NAME][cat] += StructuralIdentifiability._runtime_logger[cat]
        end
        # for cat in GB_TIME_CATEGORIES
        #     runtimes[NAME][cat] += ParamPunPam._runtime_logger[cat]
        # end
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

    open((@__DIR__) * "/systems/$NAME/timings", "w") do io
        write(io, timings)
    end
end

function dump_results()
    open((@__DIR__) * "/systems/$NAME/id_funcs", "w") do io
        if haskey(results, NAME)
            println(io, results[NAME])
        end
    end
end

process_system()
dump_timings()
dump_results()

# function generate_system()
#     idx = findfirst(b -> b[:name] == NAME, benchmarks)
#     bmark = benchmarks[idx]
#     name = bmark[:name]
#     @info "Processing $name"
#     runtimes[bmark[:name]] = Dict(c => 0.0 for c in ALL_CATEGORIES)
#     mkpath((@__DIR__)*"/systems/$name/")
#     fd = open((@__DIR__)*"/systems/$name/$name.jl", "w")
#     println(fd, "# $name")
#     println(fd, "#! format: off")
#     println(fd, "using StructuralIdentifiability")
#     println(fd, "")
#     ode = join(map(s -> "\t"*s, split(repr(bmark[:ode]), "\n")[1:end-1]), ",\n")
#     println(fd, "system = @ODEmodel(\n$ode\n)")
#     close(fd)
# end
