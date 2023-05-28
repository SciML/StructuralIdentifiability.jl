using StructuralIdentifiability

using StructuralIdentifiability.Test
using StructuralIdentifiability.TestSetExtensions

using StructuralIdentifiability.Nemo
using StructuralIdentifiability.ModelingToolkit
using StructuralIdentifiability:
    check_field_membership,
    check_identifiability,
    check_primality_zerodim,
    check_primality,
    det_minor_expansion,
    ExpVectTrie,
    get_max_below,
    ps_ode_solution,
    power_series_solution,
    ps_diff,
    ps_integrate,
    ps_matrix_inv,
    ps_matrix_homlinear_de,
    ps_matrix_linear_de,
    ps_matrix_log,
    reduce_ode_mod_p,
    ODE,
    @ODEmodel,
    truncate_matrix,
    find_ioequations,
    str_to_var,
    unpack_fraction,
    assess_global_identifiability,
    differentiate_output,
    var_to_str,
    switch_ring,
    eval_at_dict,
    assess_local_identifiability,
    assess_identifiability,
    monomial_compress,
    parent_ring_change,
    ps_matrix_const_term,
    decompose_derivative,
    PBRepresentation,
    find_leader,
    common_ring,
    lc_univariate,
    pseudodivision,
    diffreduce,
    io_switch!,
    add_outputs,
    find_ioprojections,
    choose,
    sequence_solution,
    differentiate_sequence_solution,
    differentiate_sequence_output,
    _assess_local_identifiability_discrete

function random_ps(ps_ring, range = 1000)
    result = zero(ps_ring)
    t = gen(ps_ring)
    for i in 0:(max_precision(ps_ring) - 1)
        result += (rand(Int) % range) * t^i
    end
    return result
end

function random_ps_matrix(ps_ring, matrix_space)
    result = zero(matrix_space)
    while StructuralIdentifiability.LinearAlgebra.det(ps_matrix_const_term(result)) == 0
        result = map(e -> random_ps(ps_ring), zero(matrix_space))
    end
    return result
end

@info "Testing started"

@testset "All the tests" begin
    @includetests ARGS
end
