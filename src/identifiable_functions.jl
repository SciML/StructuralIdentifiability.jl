
mutable struct IdentifiablePool
    initial_funcs::Any
    simple_funcs::Any
    simplest_funcs::Any
    all_funcs::Any
    meta::Any
    runtimes::Any

    function IdentifiablePool(id_funcs::Vector{T}, fan::FieldGeneratorsFan) where {T}
        runtimes = copy(_runtime_logger)
        meta = copy(fan.meta)
        newmeta = []
        for i in 1:length(meta)
            push!(newmeta, (meta[i].runtime, ordering = fan.orderings[i]))
        end
        new(
            id_funcs,
            fan.generators,
            fan.simplest_generators,
            fan.all_funcs,
            fan.meta,
            runtimes,
        )
    end
end

function simple_generating_set(pool::IdentifiablePool)
    pool.simplest_funcs
end

function Base.show(io::IO, pool::IdentifiablePool)
    initial_funcs = pool.initial_funcs
    funcs = pool.simple_funcs
    simplest_funcs = pool.simplest_funcs
    all_funcs = pool.all_funcs
    meta = pool.meta
    n = length(funcs)
    println(io, "Pool of identifiabile functions")
    println(io, "Initial generators: $(sum(length, initial_funcs)) functions")
    println(
        io,
        "Computed $(n) additional sets of generators in $(pool.runtimes[:id_total]) s",
    )
    println(io, "============================================")
    println(io, "The simplest set: $(length(simplest_funcs)) functions")
    println(io, "[")
    for f in simplest_funcs
        println(io, "\t", f, ",")
    end
    println(io, "]")
    println(io, "============================================")
    println(io, "All discovered functions:")
    println(io, "[")
    for f in all_funcs
        println(io, "\t", f, ",")
    end
    println(io, "]")
    println(io, "")
end

"""
    initial_identifiable_functions(ode; options...)

Returns the initial non-simplified identifiabile functions of the given ODE
system. These are presented by the coefficients of the IO-equations.

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

"""
    find_identifiable_functions(ode::ODE; options...)

Finds all functions that are identifiable for the given ODE system.

## Options

This functions takes the following optional arguments:
- `p`: A float in the range from 0 to 1, the probability of correctness. Default
  is `0.99`.
- `simplify`: When `true`, tries to simplify the identifiabile functions, and
    returns a minimal algebraically indepdendent set of functions. Default is
    `true`.
- `with_states`: Also report the identifiabile functions in states. Default is
    `false`.

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
    _runtime_logger[:id_nemo_factor] = 0.0
    _runtime_logger[:id_primality_evaluate] = 0.0
    _runtime_logger[:id_certain_factors] = []
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
        id_funcs_fracs =
            simplified_generating_set(RationalFunctionField(id_funcs), p = p, seed = seed)
    else
        id_funcs_fracs = dennums_to_fractions(id_funcs)
    end
    _runtime_logger[:id_total] = (time_ns() - runtime_start) / 1e9
    @info "The search for identifiable functions concluded in $(_runtime_logger[:id_total]) seconds"
    id_funcs_fracs
end

function describe_identifiable_pool(
    ode::ODE{T};
    p::Float64 = 0.99,
    seed = 42,
    with_states = false,
) where {T <: MPolyElem{fmpq}}
    runtime_start = time_ns()
    _runtime_logger[:id_uncertain_factorization] = 0.0
    _runtime_logger[:id_nemo_factor] = 0.0
    _runtime_logger[:id_primality_evaluate] = 0.0
    _runtime_logger[:id_certain_factors] = []
    id_funcs = initial_identifiable_functions(
        ode,
        with_states = with_states,
        adjoin_identifiable = true,
    )
    fan = generating_sets_fan(RationalFunctionField(id_funcs), p = p, seed = seed)
    _runtime_logger[:id_total] = (time_ns() - runtime_start) / 1e9
    id_pool = IdentifiablePool(id_funcs, fan)
    id_pool
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
