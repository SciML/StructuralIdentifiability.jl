# Logging used by StructuralIdentifiability.jl
# Consists of two parts:
#   - si_logger (using Logging package) for logs
#   - runtime_logger dictionary for recording runtime info
#
#   Defines functions to refresh both + specific logging functions
#   info_si, debug_si, warn_si, and error_si calling @info, @debug, @warn, and @error, respectively

# not the exhaustive list, just the ones to be initialized
const _runtime_rubrics = (
    :id_calls_to_gb,
    :id_groebner_time,
    :id_inclusion_check_mod_p,
    :id_npoints_degree,
    :id_npoints_interpolation,
    :id_beautifulization,
)

const _runtime_logger = Dict(
    :id_calls_to_gb => 0,
    :id_groebner_time => 0.0,
    :id_inclusion_check_mod_p => 0,
    :id_npoints_degree => 0,
    :id_npoints_interpolation => 0,
    :id_beautifulization => 0,
)

const _si_logger =
    Ref{Logging.ConsoleLogger}(Logging.ConsoleLogger(Logging.Info, show_limited = false))

function restart_logging(; loglevel = Logging.Info)
    _si_logger[] = Logging.ConsoleLogger(loglevel, show_limited = false)
    for r in _runtime_rubrics
        _runtime_logger[r] = 0
    end
end

function info_si(s)
    with_logger(_si_logger[]) do
        @info s
        flush(stdout)
    end
end

function debug_si(s)
    with_logger(_si_logger[]) do
        @debug s
        flush(stdout)
    end
end

function warn_si(s)
    with_logger(_si_logger[]) do
        @warn s
        flush(stdout)
    end
end

function error_si(s)
    with_logger(_si_logger[]) do
        @error s
        flush(stdout)
    end
end
