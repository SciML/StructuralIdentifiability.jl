"""
    extract_identifiable_parameters(io_equations, ode)

    Takes as input input-output equations and the corresponding ode.
    Performs the Wronskian check and extract the coefficients.

    Returns a list of the coefficients of the io_equations
"""
@timeit _to function extract_identifiable_parameters(
        io_equations::Dict{P, P},
        ode::ODE{P},
    ) where {P <: MPolyRingElem{QQFieldElem}}

    result = []

    if isempty(ode.parameters)
        return result
    end

    @info "Computing Wronskians"
    wrnsk_time = @elapsed wrnsk = wronskian(io_equations, ode)
    @info "Computed Wronskians in $wrnsk_time seconds"
    _runtime_logger[:wrnsk_time] = wrnsk_time

    dims = map(ncols, wrnsk)
    @info "Dimensions of the Wronskians $dims"

    rank_times = @elapsed wranks = map(rank, wrnsk)
    @debug "Dimensions of the Wronskians $dims"
    @debug "Ranks of the Wronskians $wranks"
    @info "Ranks of the Wronskians computed in $rank_times seconds"
    _runtime_logger[:rank_time] = rank_times

    if any([dim != rk + 1 for (dim, rk) in zip(dims, wranks)])
        @warn "One of the Wronskians has corank greater than one, so the results of the algorithm will be valid only for multiexperiment identifiability. If you still would like to assess single-experiment identifiability, we recommend using SIAN (https://github.com/alexeyovchinnikov/SIAN-Julia) or transforming all the parameters to states with zero derivative"
    end

    @debug "Extracting coefficients"
    nonparameters = filter(
        v -> !(var_to_str(v) in map(var_to_str, ode.parameters)),
        gens(parent(first(values(io_equations)))),
    )
    for eq in values(io_equations)
        eq_coefs = collect(values(extract_coefficients(eq, nonparameters)))
        sort!(eq_coefs, by = (p -> (total_degree(p), length(p))))
        append!(result, [c // eq_coefs[1] for c in eq_coefs[2:end]])
    end

    return result
end

# ------------------------------------------------------------------------------

"""
    initial_identifiable_functions(ode; options...)

Returns the partially simplified identifiabile functions of the given ODE system.
These are presented by the coefficients of the IO-equations.

## Options

This function takes the following optional arguments:
- `prob_threshold`: probability of correctness
- `with_states`: Also report the identifiabile functions in states. Default is
  `maybe` meaning that the states may be added if necessary. 
  Other options are `yes` and `no`.
  If the states are added, the identifiable functions involving parameters only
  will be simplified
- `output_saturation` - if the output saturation should be applied to models with inputs (`true` by default)

The function returns a list of identifiable functions
"""
@timeit _to function initial_identifiable_functions(
        ode::ODE{T};
        prob_threshold::Float64,
        with_states = :maybe,
        var_change_policy = :default,
        rational_interpolator = :VanDerHoevenLecerf,
        cmp = cmp_default(ode),
        output_saturation = true,
    ) where {T}
    @assert with_states in (:no, :yes, :maybe)
    ode = preprocess_ode(ode, output_saturation = output_saturation)

    @info "Computing IO-equations"
    ioeq_time = @elapsed io_equations, is_prime =
        _find_ioequations(ode; var_change_policy = var_change_policy)
    if !is_prime && with_states == :no
        @warn "Identifiable functions with parameters only requested explicitly, computing an extra projection. Note that computing the whole observation field (i.e. identifiable functions in parameters and states) may be faster"
        @info "Computing extra projection"
        io_equations = find_extra_projection(ode, io_equations)
    end
    @debug "Sizes: $(map(length, values(io_equations)))"
    @info "Computed IO-equations in $ioeq_time seconds"
    _runtime_logger[:id_io_time] = ioeq_time

    id_funcs = extract_identifiable_parameters(io_equations, ode)

    add_states = (with_states == :yes) || (with_states == :maybe && !is_prime)

    if add_states && !isempty(ode.parameters)
        @debug "Generators of identifiable functions involve states or the ideal of IO-projections is not prime, the parameter-only part is getting simplified"
        _runtime_logger[:check_time] =
            @elapsed id_funcs = simplified_generating_set(
            RationalFunctionField(id_funcs),
            prob_threshold = prob_threshold,
            seed = 42,
            simplify = :standard,
            rational_interpolator = rational_interpolator,
            cmp = cmp,
        )
    end

    if add_states
        @debug "Computing Lie derivatives"
        append!(id_funcs, states_generators(ode, io_equations))
    end

    bring, _ = Nemo.polynomial_ring(
        base_ring(ode.poly_ring),
        [var_to_str(v) for v in vcat(ode.x_vars, ode.parameters)]
    )
    id_funcs = [parent_ring_change(f, bring) for f in id_funcs]

    return id_funcs
end

# ------------------------------------------------------------------------------
"""
    check_identifiability(ode, funcs_to_check; known, prob_threshold, var_change_policy)

Input:
- `ode` - the ODE model
- `funcs_to_check` - the functions to check identifiability for
- `known` - a list of functions in states which are assumed to be known and generic
- `prob_threshold` - probability of correctness
- `var_change` - a policy for variable change (`:default`, `:yes`, `:no`), affects only the runtime

Output: a list L of booleans with L[i] being the identifiability status of the i-th function to check
"""
@timeit _to function check_identifiability(
        ode::ODE{P},
        funcs_to_check::Array{<:Any, 1};
        known::Array{P, 1} = Array{P, 1}(),
        prob_threshold::Float64 = 0.99,
        var_change_policy = :default,
    ) where {P <: MPolyRingElem{QQFieldElem}}

    states_needed = any([!is_parametric(f, ode) for f in vcat(funcs_to_check, known)])
    states_needed && @debug "Functions to check or known functions involve states"

    half_p = 0.5 + prob_threshold / 2
    identifiable_functions = initial_identifiable_functions(
        ode,
        prob_threshold = half_p,
        var_change_policy = var_change_policy,
        with_states = states_needed ? :yes : :maybe,
    )

    if isempty(identifiable_functions)
        return [isempty(vars(f)) for f in funcs_to_check]
    end

    bring = parent(numerator(first(identifiable_functions)))
    append!(identifiable_functions, [parent_ring_change(f, bring) // one(bring) for f in known])
    funcs_to_check = Vector{Generic.FracFieldElem{P}}(
        map(f -> parent_ring_change(f, bring) // one(bring), funcs_to_check),
    )

    _runtime_logger[:check_time] =
        get(_runtime_logger, :check_time, 0.0) + @elapsed result = field_contains(
        RationalFunctionField(identifiable_functions),
        funcs_to_check,
        half_p,
    )
    return result
end

function check_identifiability(
        ode::ODE{P};
        known::Array{P, 1} = Array{P, 1}(),
        prob_threshold::Float64 = 0.99,
        var_change_policy = :default,
    ) where {P <: MPolyRingElem{QQFieldElem}}
    return check_identifiability(
        ode,
        ode.parameters,
        known = known,
        prob_threshold = prob_threshold,
        var_change_policy = var_change_policy,
    )
end

#------------------------------------------------------------------------------
"""
    assess_global_identifiability(ode::ODE{P}, prob_threshold::Float64=0.99; var_change=:default) where P <: MPolyRingElem{QQFieldElem}

Input:
- `ode` - the ODE model
- `known` - a list of functions in states which are assumed to be known and generic
- `prob_threshold` - probability of correctness
- `var_change` - a policy for variable change (`:default`, `:yes`, `:no`), affects only the runtime

Output:
- a dictionary mapping each parameter to a boolean.

Checks global identifiability for parameters of the model provided in `ode`. Call this function to check global identifiability of all parameters automatically.
"""
function assess_global_identifiability(
        ode::ODE{P},
        known::Array{P, 1} = Array{P, 1}(),
        prob_threshold::Float64 = 0.99;
        var_change = :default,
    ) where {P <: MPolyRingElem{QQFieldElem}}
    result_list = assess_global_identifiability(
        ode,
        ode.parameters,
        known,
        prob_threshold;
        var_change = var_change,
    )

    return Dict(param => val for (param, val) in zip(ode.parameters, result_list))
end

#------------------------------------------------------------------------------

"""
    assess_global_identifiability(ode, [funcs_to_check, prob_threshold=0.99, var_change=:default])

Input:
- `ode` - the ODE model
- `funcs_to_check` - rational functions in parameters
- `known` - function in parameters that are assumed to be known and generic
- `prob_threshold` - probability of correctness
- `var_change` - a policy for variable change (`:default`, `:yes`, `:no`),
                affects only the runtime

Output:
- array of length `length(funcs_to_check)` with true/false values for global identifiability
        or dictionary `param => Bool` if `funcs_to_check` are not given

Checks global identifiability of functions of parameters specified in `funcs_to_check`.
"""
@timeit _to function assess_global_identifiability(
        ode::ODE{P},
        funcs_to_check::Array{<:Any, 1},
        known::Array{P, 1} = Array{P, 1}(),
        prob_threshold::Float64 = 0.99;
        var_change = :default,
    ) where {P <: MPolyRingElem{QQFieldElem}}
    result = check_identifiability(
        ode,
        funcs_to_check,
        known = known,
        prob_threshold = prob_threshold,
        var_change_policy = var_change,
    )

    return result
end
