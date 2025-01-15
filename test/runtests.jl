using StructuralIdentifiability

using Test
using TestSetExtensions
using SpecialFunctions

using StructuralIdentifiability.DataStructures
using StructuralIdentifiability.Nemo
using StructuralIdentifiability:
    field_contains,
    check_identifiability,
    check_primality_zerodim,
    check_primality,
    check_algebraicity,
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
    _assess_local_identifiability_discrete_aux,
    extract_coefficients,
    extract_coefficients_ratfunc,
    lie_derivative,
    states_generators,
    RationalFunctionField,
    replace_with_ic,
    x_vars,
    y_vars,
    x_equations,
    y_equations,
    inputs,
    quotient_basis,
    rational_function_cmp

const GROUP = get(ENV, "GROUP", "All")

@static if VERSION >= v"1.10.0"
    if GROUP == "All" || GROUP == "ModelingToolkitSIExt"
        using Pkg
        Pkg.add("ModelingToolkit")
        Pkg.add("Symbolics")
        using ModelingToolkit
        using Symbolics
    end
end

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

function get_test_files(group)
    result = Vector{String}()
    for (dir, _, files) in walkdir("./")
        for fname in files
            if fname != "runtests.jl" && endswith(fname, ".jl")
                if group == "All" ||
                   (group == "Core" && dir != "./extensions") ||
                   (
                       group == "ModelingToolkitSIExt" &&
                       dir == "./extensions" &&
                       VERSION >= v"1.10.0"
                   )
                    push!(result, dir * "/" * fname)
                end
            end
        end
    end
    return result
end

@info "Testing started"

all_tests = get_test_files(GROUP)
if !isempty(ARGS)
    all_tests = ARGS
end

@time @testset "All the tests" verbose = true begin
    for test_file in all_tests
        include(test_file)
    end
end
