using StructuralIdentifiability

using StructuralIdentifiability.Test
using StructuralIdentifiability.TestSetExtensions

using StructuralIdentifiability.AbstractAlgebra
using StructuralIdentifiability.Nemo
using StructuralIdentifiability.Singular

using StructuralIdentifiability: check_field_membership, check_identifiability, check_primality_zerodim,
                                 det_minor_expansion, ExpVectTrie, get_max_below, ps_ode_solution, 
                                 power_series_solution, ps_diff, ps_integrate, ps_matrix_inv,
                                 ps_matrix_homlinear_de, ps_matrix_linear_de, ps_matrix_log,
                                 reduce_ode_mod_p, simplify_field_generators, ODE, @ODEmodel,
                                 truncate_matrix, find_ioequations, str_to_var, unpack_fraction,
                                 assess_global_identifiability, differentiate_output, var_to_str,
                                 switch_ring, eval_at_dict, assess_local_identifiability,
                                 assess_identifiability, monomial_compress, parent_ring_change

function random_ps(ps_ring, range = 1000)
    result = zero(ps_ring)
    t = gen(ps_ring)
    for i in 0:(max_precision(ps_ring) - 1)
        result += (rand(Int) % range) * t^i
    end
    return result
end

function random_ps_matrix(ps_ring, matrix_space)
    return map(e -> random_ps(ps_ring), zero(matrix_space))
end


@info "Testing started"

@testset "All the tests" begin
    @includetests ARGS
end

