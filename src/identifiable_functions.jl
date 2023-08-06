
"""
    find_identifiable_functions(ode::ODE; options...)

Finds all functions of parameters/states that are identifiable in the given ODE
system.

## Options

This functions takes the following optional arguments:
- `p`: A float in the range from 0 to 1, the probability of correctness. Default
  is `0.99`.
- `simplify`: When `true`, tries to simplify the identifiabile functions, and
    returns a minimal algebraically indepdendent set. Default is `true`.
- `with_states`: When `true`, also reports the identifiabile functions in the
    ODE states. Default is `false`.
- `strategy`: The simplification strategy. Possible options are:
    - `(:gb, )`: Extract the coefficients of a Groebner basis of the MQS ideal
      (default).
    - `(:gbfan, N)`: Same as `:gb`, but computes `N` bases for different
        rankings of variables.
    - `(:normalforms, N)`: Same as `:gb`, but adjoins the results of normal form
      computations of monomials up to the total degree `N`. 
    - `(:hybrid, )`: The best of all worlds.
- `seed`: The rng seed. Default value is `42`.

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
    adjoin_identifiable = false,
    strategy = (:gb,),
) where {T <: MPolyElem{fmpq}}
    @assert first(strategy) in (:gb, :gbfan, :normalforms, :hybrid)
    runtime_start = time_ns()
    id_funcs = initial_identifiable_functions(
        ode,
        with_states = with_states,
        adjoin_identifiable = adjoin_identifiable,
    )
    if simplify
        if isempty(id_funcs)
            bring = parent(ode)
            id_funcs = [one(bring)]
        end
        id_funcs_fracs = simplified_generating_set(
            RationalFunctionField(id_funcs),
            p = p,
            seed = seed,
            strategy = strategy,
        )
    else
        id_funcs_fracs = dennums_to_fractions(id_funcs)
    end
    _runtime_logger[:id_total] = (time_ns() - runtime_start) / 1e9
    @info "The search for identifiable functions concluded in $(_runtime_logger[:id_total]) seconds"
    id_funcs_fracs
end

"""
    find_identifiable_functions(ode::ModelingToolkit.ODESystem; measured_quantities=[], options...)

Finds all functions that are identifiable for the given ODE system.
    
## Options

This functions takes the following optional arguments:
- `measured_quantities` - the output functions of the model.

## Example

```jldoctest
using StructuralIdentifiability
using ModelingToolkit

@parameters a01 a21 a12
@variables t x0(t) x1(t) y1(t) [output = true]
D = Differential(t)

eqs = [
    D(x0) ~ -(a01 + a21) * x0 + a12 * x1, 
    D(x1) ~ a21 * x0 - a12 * x1, y1 ~ x0
]
de = ODESystem(eqs, t, name = :Test)

find_identifiable_functions(de, measured_quantities = [y1 ~ x0])

# prints
2-element Vector{Num}:
         a01*a12
 a01 + a12 + a21
```
"""
function find_identifiable_functions(
    ode::ModelingToolkit.ODESystem;
    measured_quantities = Array{ModelingToolkit.Equation}[],
    simplify = true,
    with_states = false,
    p = 0.99,
    seed = 42,
    strategy = (:gb,),
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
        strategy = strategy,
    )
    nemo2mtk = Dict(v => Num(k) for (k, v) in conversion)
    out_funcs = [eval_at_dict(func, nemo2mtk) for func in result]
    return out_funcs
end

"""
    initial_identifiable_functions(ode; options...)

Returns the non-simplified identifiabile functions of the given ODE system.
These are presented by the coefficients of the IO-equations.

## Options

This functions takes the following optional arguments:
- `with_states`: Also report the identifiabile functions in states. Default is
  `false`.
- `adjoin_identifiable`: Adjoin globally identifiabile parameters to the output.
  *Global identifiability is assessed modulo a prime.* Default is `true`.
"""
function initial_identifiable_functions(
    ode::ODE{T};
    with_states = false,
    adjoin_identifiable = true,
) where {T}
    @info "Computing IO-equations"
    runtime = @elapsed io_equations = find_ioequations(ode)
    @info "IO-equations computed in $runtime seconds"
    _runtime_logger[:id_io_time] = runtime
    id_funcs = extract_identifiable_functions_raw(
        io_equations,
        ode,
        empty(ode.parameters),
        with_states,
    )
    _runtime_logger[:id_global_time] = 0.0
    if adjoin_identifiable
        @info "Assessing global identifiability"
        to_check = ode.parameters
        if with_states
            to_check = vcat(to_check, ode.x_vars)
        end
        bring = parent(first(first(id_funcs)))
        to_check = map(f -> parent_ring_change(f, bring), to_check)
        runtime = @elapsed global_result = check_field_membership_mod_p(id_funcs, to_check)
        @debug "Identifiable parameters and states are" to_check global_result
        @info "Global identifiability assessed in $runtime seconds"
        _runtime_logger[:id_global_time] = runtime
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
        id_funcs = vcat(id_funcs, known_quantities)
    end
    return id_funcs
end
