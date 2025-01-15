
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
- `known_ic`: a list of functions whose initial conditions are assumed to be known,
  then the returned identifiable functions will be functions of parameters and
  initial conditions, not states (this is an experimental functionality).
- `prob_threshold`: A float in the range from 0 to 1, the probability of correctness. Default
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
3-element Vector{AbstractAlgebra.Generic.FracFieldElem{Nemo.QQMPolyRingElem}}:
 a12 + a01 + a21
 a12*a01
```

"""
function find_identifiable_functions(
    ode::ODE{T};
    known_ic::Vector{<:ExtendedFraction{T}} = Vector{ExtendedFraction{T}}(),
    prob_threshold::Float64 = 0.99,
    seed = 42,
    with_states = false,
    simplify = :standard,
    rational_interpolator = :VanDerHoevenLecerf,
    loglevel = Logging.Info,
) where {T <: MPolyRingElem{QQFieldElem}}
    restart_logging(loglevel = loglevel)
    reset_timings()
    with_logger(_si_logger[]) do
        if isempty(known_ic)
            return _find_identifiable_functions(
                ode,
                prob_threshold = prob_threshold,
                seed = seed,
                with_states = with_states,
                simplify = simplify,
                rational_interpolator = rational_interpolator,
            )
        else
            id_funcs = _find_identifiable_functions_kic(
                ode,
                known_ic,
                prob_threshold = prob_threshold,
                seed = seed,
                simplify = simplify,
                rational_interpolator = rational_interpolator,
            )
            # renaming variables from `x(t)` to `x(0)`
            return replace_with_ic(ode, id_funcs)
        end
    end
end

function _find_identifiable_functions(
    ode::ODE{T};
    prob_threshold::Float64 = 0.99,
    seed = 42,
    with_states = false,
    simplify = :standard,
    rational_interpolator = :VanDerHoevenLecerf,
) where {T <: MPolyRingElem{QQFieldElem}}
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
    half_p = 0.5 + prob_threshold / 2
    id_funcs, bring = initial_identifiable_functions(
        ode,
        prob_threshold = half_p,
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
            prob_threshold = half_p,
            seed = seed,
            simplify = simplify,
            rational_interpolator = rational_interpolator,
            priority_variables = [parent_ring_change(p, bring) for p in ode.parameters],
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
