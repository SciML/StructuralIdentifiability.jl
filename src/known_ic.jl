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
    known_ic::Vector{<:Union{T, Generic.Frac{T}}};
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
    known_ic::Vector{<:Union{T, Generic.Frac{T}}};
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

    id_funcs = simplified_generating_set(
        RationalFunctionField(
            vcat(id_funcs_general, [f // one(parent(ode)) for f in known_ic]),
        ),
        p = half_p,
        seed = seed,
        simplify = simplify,
        rational_interpolator = rational_interpolator,
    )

    @info "The search for identifiable functions with known initial conditions concluded in $((time_ns() - runtime_start) / 1e9) seconds"

    return replace_with_ic(ode, id_funcs)
end

"""
    assess_identifiability_kic(ode, known_ic; funcs_to_check = [], p=0.99, loglevel=Logging.Info)

Input:
- `ode` - the ODE model
- `known_ic` - a list of functions for which initial conditions are assumed to be known and generic
- `funcs_to_check` - list of functions to check identifiability for; if empty, all parameters
   and states are taken
- `p` - probability of correctness.
- `loglevel` - the minimal level of log messages to display (`Logging.Info` by default)

Assesses identifiability of parameters and initial conditions of a given ODE model. 
The result is guaranteed to be correct with the probability at least `p`.
The function returns an (ordered) dictionary from the functions to check to their identifiability properties 
(one of `:nonidentifiable`, `:locally`, `:globally`).
"""
function assess_identifiability_kic(
    ode::ODE{P},
    known_ic::Vector{<:Union{P, Generic.Frac{P}}};
    funcs_to_check = Vector(),
    p::Float64 = 0.99,
    loglevel = Logging.Info,
) where {P <: MPolyElem{fmpq}}
    restart_logging(loglevel = loglevel)
    reset_timings()
    with_logger(_si_logger[]) do
        return _assess_identifiability_kic(
            ode,
            known_ic,
            p = p,
            funcs_to_check = funcs_to_check,
        )
    end
end

function _assess_identifiability_kic(
    ode::ODE{P},
    known_ic::Vector{<:Union{P, Generic.Frac{P}}};
    funcs_to_check = Vector(),
    p::Float64 = 0.99,
) where {P <: MPolyElem{fmpq}}
    runtime_start = time_ns()
    if length(funcs_to_check) == 0
        funcs_to_check = vcat(ode.x_vars, ode.parameters)
    end
    half_p = 0.5 + p / 2
    id_funcs = _find_identifiable_functions_kic(ode, known_ic, p = half_p)
    funcs_to_check = replace_with_ic(ode, funcs_to_check)
    result = OrderedDict(f => :globally for f in funcs_to_check)

    half_p = 0.5 + p / 2
    local_result = check_algebraicity(
        RationalFunctionField([[denominator(f), numerator(f)] for f in id_funcs]),
        [f // one(parent(f)) for f in funcs_to_check],
        half_p,
    )
    global_result = field_contains(
        RationalFunctionField([[denominator(f), numerator(f)] for f in id_funcs]),
        [f // one(parent(f)) for f in funcs_to_check],
        half_p,
    )
    for (i, f) in enumerate(funcs_to_check)
        if !local_result[i]
            result[f] = :nonidentifiable
        elseif !global_result[i]
            result[f] = :locally
        end
    end
    @info "Assessing identifiability with known initial conditions concluded in $((time_ns() - runtime_start) / 1e9) seconds"
    return result
end
