
"""
    find_identifiable_functions(ode::ODE; options...)

Finds all functions of parameters/states that are identifiable in the given ODE
system.

## Options

This functions takes the following optional arguments:
- `with_states`: When `true`, also reports the identifiabile functions in the
    ODE states. Default is `false`.
- `simplify`: The extent to which the output functions are simplified. Stronger
  simplification may require more time. Possible options are:
  - `:standard`: Default simplification.
  - `:weak`: Weak simplification. This option is the fastest, but the output
    functions can be quite complex.
  - `:strong`: Strong simplification. This option is the slowest, but the output
  functions are nice and simple.
  - `:absent`: No simplification.
- `p`: A float in the range from 0 to 1, the probability of correctness. Default
  is `0.99`.
- `seed`: The rng seed. Default value is `42`.
- `loglevel` - the minimal level of log messages to display (`Logging.Info` by default)

## Example

```jldoctest
using StructuralIdentifiability

ode = @ODEmodel(
    x0'(t) = -(a01 + a21) * x0(t) + a12 * x1(t),
    x1'(t) = a21 * x0(t) - a12 * x1(t),
    y(t) = x0(t)
)

find_identifiable_functions(ode)

# prints
3-element Vector{AbstractAlgebra.Generic.Frac{Nemo.fmpq_mpoly}}:
 a12 + a01 + a21
 a12*a01
```

"""
function find_identifiable_functions(
    ode::ODE{T};
    p::Float64 = 0.99,
    seed = 42,
    with_states = false,
    simplify = :standard,
    rational_interpolator = :VanDerHoevenLecerf,
    loglevel = Logging.Info,
) where {T <: MPolyElem{fmpq}}
    restart_logging(loglevel = loglevel)
    reset_timings()
    with_logger(_si_logger[]) do
        return _find_identifiable_functions(
            ode,
            p = p,
            seed = seed,
            with_states = with_states,
            simplify = simplify,
            rational_interpolator = rational_interpolator,
        )
    end
end

function _find_identifiable_functions(
    ode::ODE{T};
    p::Float64 = 0.99,
    seed = 42,
    with_states = false,
    simplify = :standard,
    rational_interpolator = :VanDerHoevenLecerf,
) where {T <: MPolyElem{fmpq}}
    Random.seed!(seed)
    @assert simplify in (:standard, :weak, :strong, :absent)
    runtime_start = time_ns()
    if isempty(ode.parameters) && !with_states
        @warn """
        There are no parameters in the given ODE, thus no identifiabile
        functions.
        Use `find_identifiable_functions` with keyword `with_states=true` to
        compute the functions with the ODE states included."""
        bring = parent(ode)
        id_funcs = [one(bring)]
        return id_funcs
    end
    half_p = 0.5 + p / 2
    id_funcs, bring = initial_identifiable_functions(
        ode,
        p = half_p,
        with_states = with_states,
        rational_interpolator = rational_interpolator,
    )
    # If simplification is needed
    if simplify !== :absent
        if isempty(id_funcs)
            bring = parent(ode)
            id_funcs = [one(bring)]
        end
        id_funcs_fracs = simplified_generating_set(
            RationalFunctionField(id_funcs),
            p = half_p,
            seed = seed,
            simplify = simplify,
            rational_interpolator = rational_interpolator,
        )
    else
        id_funcs_fracs = dennums_to_fractions(id_funcs)
    end
    id_funcs_fracs = [parent_ring_change(f, parent(ode)) for f in id_funcs_fracs]
    _runtime_logger[:id_total] = (time_ns() - runtime_start) / 1e9
    _runtime_logger[:are_id_funcs_polynomial] = all(isone ∘ denominator, id_funcs_fracs)
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
- `loglevel` - the verbosity of the logging 
  (can be Logging.Error, Logging.Warn, Logging.Info, Logging.Debug)

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
    p::Float64 = 0.99,
    seed = 42,
    with_states = false,
    simplify = :standard,
    rational_interpolator = :VanDerHoevenLecerf,
    loglevel = Logging.Info,
)
    restart_logging(loglevel = loglevel)
    reset_timings()
    with_logger(_si_logger[]) do
        return _find_identifiable_functions(
            ode,
            measured_quantities = measured_quantities,
            p = p,
            seed = seed,
            with_states = with_states,
            simplify = simplify,
            rational_interpolator = rational_interpolator,
        )
    end
end

function _find_identifiable_functions(
    ode::ModelingToolkit.ODESystem;
    measured_quantities = Array{ModelingToolkit.Equation}[],
    p::Float64 = 0.99,
    seed = 42,
    with_states = false,
    simplify = :standard,
    rational_interpolator = :VanDerHoevenLecerf,
)
    Random.seed!(seed)
    if isempty(measured_quantities)
        measured_quantities = get_measured_quantities(ode)
    end
    ode, conversion = mtk_to_si(ode, measured_quantities)
    result = _find_identifiable_functions(
        ode,
        simplify = simplify,
        p = p,
        seed = seed,
        with_states = with_states,
        rational_interpolator = rational_interpolator,
    )
    result = [parent_ring_change(f, ode.poly_ring) for f in result]
    nemo2mtk = Dict(v => Num(k) for (k, v) in conversion)
    out_funcs = [eval_at_dict(func, nemo2mtk) for func in result]
    return out_funcs
end
