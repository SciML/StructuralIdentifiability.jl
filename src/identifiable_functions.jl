
"""
    find_identifiable_functions(ode::ODE; options...)

Finds all functions of parameters/states that are identifiable in the given ODE
system.

## Options

This functions takes the following optional arguments:
- `p`: A float in the range from 0 to 1, the probability of correctness. Default
  is `0.99`.
- `simplify`: When `true`, tries to simplify the identifiabile functions, and
    returns an algebraically indepdendent set. Default is `true`.
- `with_states`: When `true`, also reports the identifiabile functions in the
    ODE states. Default is `false`.
- `strategy`: The simplification strategy. Possible options are:
    - `(:gb, )`: Extract the coefficients of a Groebner basis of the MQS ideal
      (default).
    - `(:gbfan, N)`: Same as `:gb`, but computes `N` bases for different random
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
    half_p = 0.5 + p / 2
    id_funcs, bring =
        initial_identifiable_functions(ode, p = half_p, with_states = with_states)
    if simplify
        if isempty(id_funcs)
            bring = parent(ode)
            id_funcs = [one(bring)]
        end
        id_funcs_fracs = simplified_generating_set(
            RationalFunctionField(id_funcs),
            p = half_p,
            seed = seed,
            strategy = strategy,
        )
    else
        id_funcs_fracs = dennums_to_fractions(id_funcs)
    end
    id_funcs_fracs = [parent_ring_change(f, parent(ode)) for f in id_funcs_fracs]
    _runtime_logger[:id_total] = (time_ns() - runtime_start) / 1e9
    @info "The search for identifiable functions concluded in $(_runtime_logger[:id_total]) seconds"
    return id_funcs_fracs
end

"""
    find_identifiable_functions(ode::ModelingToolkit.ODESystem; measured_quantities=[], options...)

Finds all functions of parameters/states that are identifiable in the given ODE
system.

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
    result = [parent_ring_change(f, ode.poly_ring) for f in result]
    nemo2mtk = Dict(v => Num(k) for (k, v) in conversion)
    out_funcs = [eval_at_dict(func, nemo2mtk) for func in result]
    return out_funcs
end
