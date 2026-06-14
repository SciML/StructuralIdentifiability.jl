using StructuralIdentifiability

using Test
using TestSetExtensions
using SpecialFunctions

using StructuralIdentifiability.DataStructures
using StructuralIdentifiability.Nemo
using StructuralIdentifiability.RationalFunctionFields
using StructuralIdentifiability.RationalFunctionFields:
    str_to_var, unpack_fraction, var_to_str, eval_at_dict, fractionfree_generators_raw
using StructuralIdentifiability:
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
    assess_global_identifiability,
    differentiate_output,
    switch_ring,
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
    _assess_local_identifiability_discrete_aux,
    extract_coefficients,
    extract_coefficients_ratfunc,
    lie_derivative,
    states_generators,
    replace_with_ic,
    x_vars,
    y_vars,
    x_equations,
    y_equations,
    inputs,
    quotient_basis,
    propose_orders,
    saturate_outputs

const GROUP = get(ENV, "GROUP", "All")

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

function rand_poly(deg, vars)
    result = 0
    indices = vcat(collect(1:length(vars)), collect(1:length(vars)))
    monomials = []
    for d in 0:deg
        for subs in StructuralIdentifiability.IterTools.subsets(indices, d)
            push!(monomials, subs)
        end
    end

    for subs in monomials
        monom = rand(-50:50)
        for v_ind in subs
            monom *= vars[v_ind]
        end
        result += monom
    end

    return result
end
