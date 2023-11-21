"""
    find_identifiable_functions_kic(ode::ODE, known_ic; options...)

Finds all functions of parameters/initial conditions that are identifiable in the given ODE
system under assumptions that the initial conditions for functions in the list
`known_ic` are known and generic.

## Options

This functions takes the following optional arguments:
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

**This is experimental functionality**

```

"""
function find_identifiable_functions_kic(
    ode::ODE{T},
    known_ic::Vector{<: Union{T, Generic.Frac{T}}};
    p::Float64 = 0.99,
    seed = 42,
    simplify = :standard,
    rational_interpolator = :VanDerHoevenLecerf,
    loglevel = Logging.Info,
) where {T <: MPolyElem{fmpq}}
    restart_logging(loglevel = loglevel)
    reset_timings()
    with_logger(_si_logger[]) do
        return _find_identifiable_functions_kic(
            ode,
            known_ic,
            p = p,
            seed = seed,
            simplify = simplify,
            rational_interpolator = rational_interpolator,
        )
    end
end

function _find_identifiable_functions_kic(
    ode::ODE{T},
    known_ic::Vector{<: Union{T, Generic.Frac{T}}};
    p::Float64 = 0.99,
    seed = 42,
    simplify = :standard,
    rational_interpolator = :VanDerHoevenLecerf,
) where {T <: MPolyElem{fmpq}}
    Random.seed!(seed)
    @assert simplify in (:standard, :weak, :strong, :absent)
    half_p = 0.5 + p / 2
    runtime_start = time_ns()
    id_funcs_general = find_identifiable_functions(
        ode,
        p = half_p,
        with_states = true,
        simplify = simplify,
        rational_interpolator = rational_interpolator,
        seed = seed,
    )
    @info "The search for identifiable functions with known initial conditions concluded in $((time_ns() - runtime_start) / 1e9) seconds"

    id_funcs = simplified_generating_set(
        RationalFunctionField(vcat(id_funcs_general, [f // one(parent(ode)) for f in known_ic])),
        p = half_p,
        seed = seed,
        simplify = simplify,
        rational_interpolator = rational_interpolator,
    )

    return id_funcs
end


