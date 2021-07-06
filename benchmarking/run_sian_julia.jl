using Logging
using SIAN

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

# to run this, one should comment the NFkB part in the benchmark file
include("benchmarks.jl")
NUM_RUNS = 5

index_to_run = 1

bmark = benchmarks[index_to_run]
@info "Processing $(bmark[:name])"
ode = bmark[:ode]

total_time = 0

for i in 1:(NUM_RUNS + 1)
    @info "Run numer $i"
    rtime = @elapsed identifiability_ode(ode, get_parameters(ode))
    @info "Runtime: $rtime"
    if i > 1
        global total_time += rtime
    end
end

@info total_time
@info "Average: $(total_time / NUM_RUNS)"
