
"""
    find_identifiable_functions(ode::ODE; p = 0.99, seed = 42, with_states = false)

Finds all expressions of parameters that are identifiable for the given ODE
system.

Input:
- `ode` - `ODE`-system.
Output:
- a set of expressions of parameters that are identifiabile

## Example

```jldoctest
using StructuralIdentifiability

de = @ODEmodel(
    x0'(t) = -(a01 + a21) * x0(t) + a12 * x1(t),
    x1'(t) = a21 * x0(t) - a12 * x1(t),
    y(t) = x0(t)
)

find_identifiable_functions(de)
# prints
3-element Vector{AbstractAlgebra.Generic.Frac{Nemo.fmpq_mpoly}}:
 a12 + a01 + a21
 a12*a01
```

"""
function find_identifiable_functions(
    ode::ODE{T};
    p::Float64 = 0.99,
    simplify = true,
    seed = 42,
    with_states = false,
    adjoin_identifiable = true,
) where {T <: MPolyElem{fmpq}}
    runtime_start = time_ns()
    _runtime_logger[:id_uncertain_factorization] = 0.0
    _runtime_logger[:id_primality_evaluate] = 0.0
    _runtime_logger[:id_certain_factors] = []

    @info "Computing IO-equations"
    runtime = @elapsed io_equations = find_ioequations(ode)
    @info "IO-equations computed in $runtime seconds"
    _runtime_logger[:id_io_time] = runtime
    identifiable_functions_raw = extract_identifiable_functions_raw(
        io_equations,
        ode,
        empty(ode.parameters),
        with_states,
    )
    id_funcs = identifiable_functions_raw
    if adjoin_identifiable
        @info "Assessing global identifiability"
        to_check = ode.parameters
        if with_states
            to_check = vcat(to_check, ode.x_vars)
        end
        bring = parent(first(first(identifiable_functions_raw)))
        to_check = map(f -> parent_ring_change(f, bring), to_check)
        runtime = @elapsed global_result =
            check_field_membership_mod_p(identifiable_functions_raw, to_check)
        @debug "Identifiable parameters and states are" to_check global_result
        @info "Global identifiability assessed in $runtime seconds"
        known_quantities = Vector{Vector{elem_type(bring)}}()
        for (glob, p) in zip(global_result, to_check)
            if glob
                if all(in.(map(var_to_str, vars(p)), [map(var_to_str, gens(bring))]))
                    push!(known_quantities, [one(bring), parent_ring_change(p, bring)])
                else
                    @warn "Known quantity $p cannot be casted and is thus dropped"
                end
            end
        end
        id_funcs = vcat(identifiable_functions_raw, known_quantities)
    end
    if simplify
        if isempty(id_funcs)
            id_funcs = [one(bring)]
        end
        id_funcs_fracs =
            simplified_generating_set(RationalFunctionField(id_funcs), p = p, seed = seed)
    else
        id_funcs_fracs = dennums_to_fractions(id_funcs)
    end
    _runtime_logger[:id_total] = (time_ns() - runtime_start) / 1e9
    @info "The search for identifiable functions concluded in $(_runtime_logger[:id_total]) seconds"
    return id_funcs_fracs
end

"""
    find_identifiable_functions(ode::ModelingToolkit.ODESystem; measured_quantities=Array{ModelingToolkit.Equation}[])

Finds all expressions of parameters that are identifiable for the given ODE
system.
    
Input:
- `ode` - the ModelingToolkit.ODESystem object that defines the model
- `measured_quantities` - the output functions of the model
Output:
- a set of expressions of parameters that are identifiabile

"""
function find_identifiable_functions(
    ode::ModelingToolkit.ODESystem;
    measured_quantities = Array{ModelingToolkit.Equation}[],
    simplify = true,
    with_states = false,
    p = 0.99,
    seed = 42,
)
    if isempty(measured_quantities)
        measured_quantities = get_measured_quantities(ode)
    end
    ode, conversion = preprocess_ode(ode, measured_quantities)
    result = find_identifiable_functions(
        ode,
        simplify = simplify,
        p = p,
        seed = seed,
        with_states = with_states,
    )
    nemo2mtk = Dict(v => Num(k) for (k, v) in conversion)
    out_funcs = [eval_at_dict(func, nemo2mtk) for func in result]
    return out_funcs
end
