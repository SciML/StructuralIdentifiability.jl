# Logging used by StructuralIdentifiability.jl
# Consists of two parts:
#   - si_logger (using Logging package) for logs
#   - runtime_logger dictionary for recording runtime info
#

# not the exhaustive list, just the ones to be initialized
const _runtime_rubrics = (
    :id_calls_to_gb,
    :id_groebner_time,
    :id_inclusion_check_mod_p,
    :id_npoints_degree,
    :id_npoints_interpolation,
    :id_beautifulization,
    :id_npoints_normalform,
    :id_normalforms_time,
)

const _runtime_logger = Dict(
    :id_calls_to_gb => 0,
    :id_groebner_time => 0.0,
    :id_inclusion_check_mod_p => 0,
    :id_npoints_degree => 0,
    :id_npoints_interpolation => 0,
    :id_beautifulization => 0,
    :id_npoints_normalform => 0,
    :id_normalforms_time => 0.0,
)

const _si_logger =
    Ref{Logging.ConsoleLogger}(Logging.ConsoleLogger(Logging.Info, show_limited = false))

function restart_logging(; loglevel = Logging.Info, stream = nothing)
    @assert loglevel isa Base.CoreLogging.LogLevel
    if stream !== nothing
    	_si_logger[] = Logging.ConsoleLogger(stream, loglevel, show_limited = false)
    else
        _si_logger[] = Logging.ConsoleLogger(loglevel, show_limited = false)
    end
    for r in _runtime_rubrics
        _runtime_logger[r] = 0
    end
    return nothing
end

###
# Timings for StructuralIdentifiability.jl.
# Uses the the @timeit macro from TimerOutputs.jl.
#
# Provides a couple of useful functions:
#   - enable_timer, for enabling or disabling the timer in SI
#   - reset_timings, for clearing the collected timings
#   - print_timings_table, for printing the table with timings

timeit_debug_enabled() = false

const _to = TimerOutputs.TimerOutput()

"""
    enable_timer(flag)

If `flag` is `true`, enables the global timer. Otherwise, disables it.
"""
function enable_timer(flag::Bool)
    if flag
        enable_timer!(_to)
    else
        disable_timer!(_to)
    end
    nothing
end

# By default, the timer is disabled
enable_timer(false)

"""
    reset_timings()

Resets the global timer.
"""
function reset_timings()
    TimerOutputs.reset_timer!(_to)
    nothing
end

"""
    print_timings_table()

Prints the table with collected timings data to `stdout`.
"""
function print_timings_table()
    iostream = stdout
    TimerOutputs.print_timer(
        iostream,
        _to,
        allocations = true,
        sortby = :time,
        linechars = :ascii,
        compact = false,
        title = "SI.jl",
    )
end

function nonrational_error(precision::String)
    return Base.ArgumentError(
        """
The system does not seem to have rational (polynomial divided by polynomial) right-hand side.
More precisely: $precision.
For the moment, such systems cannot be handled directly. We advise to try a variable transformation,
for an example and some guidelines, we refer to pages 4-5 of https://doi.org/10.3390/v17040496
(in particular the discussion after Proposition 1). Further examples of transformations can be found
in Sections A.2 and A.3 of the Supplementary Material of https://doi.org/10.1093/bioinformatics/bty1069
""",
    )
end
