using StructuralIdentifiability, BenchmarkTools

const SUITE = BenchmarkGroup()

# =============================================================================
# Model construction — linear compartment models
# =============================================================================

SUITE["construction"] = BenchmarkGroup()

graph3 = [[2], [1, 3], [1, 2]]
graph4 = [[2, 4], [1, 3], [2, 4], [1, 3]]

SUITE["construction"]["lincomp_3"] = @benchmarkable linear_compartment_model(
    $graph3; inputs = [1], outputs = [1], leaks = []
)
SUITE["construction"]["lincomp_4"] = @benchmarkable linear_compartment_model(
    $graph4; inputs = [2], outputs = [1], leaks = [2, 3]
)

# =============================================================================
# Identifiability assessment — the core workload: differential algebra +
# rational function field computation
# =============================================================================

SUITE["identifiability"] = BenchmarkGroup()

ode3 = linear_compartment_model(graph3; inputs = [1], outputs = [1])
ode4 = linear_compartment_model(graph4; inputs = [2], outputs = [1], leaks = [2, 3])

SUITE["identifiability"]["local_3"] = @benchmarkable assess_local_identifiability($ode3)
SUITE["identifiability"]["local_4"] = @benchmarkable assess_local_identifiability($ode4)
SUITE["identifiability"]["global_3"] = @benchmarkable assess_identifiability($ode3) seconds = 600

# =============================================================================
# IO-equation derivation and identifiable function extraction
# =============================================================================

SUITE["ioequations"] = BenchmarkGroup()

SUITE["ioequations"]["find_3"] = @benchmarkable find_ioequations($ode3)
SUITE["ioequations"]["find_4"] = @benchmarkable find_ioequations($ode4)
