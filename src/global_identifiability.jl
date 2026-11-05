"""
    extract_identifiable_functions_raw(io_equations, ode, known, with_states)

Takes as input input-output equations, the corresponding ode, a list of functions assumed to be known
and a flag `with_states`.
Extracts generators of the field of identifiable functions (with or without states) without
any simplifications.

Returns a tuple consisting of
- a dictionary with at most two keys: `:no_states` and `:with_states`. The respective values
are identifiable functions containing or not the state variables
- a polynomial ring containing all returned identifiable functions
(parameters or parameters + states)
"""
@timeit _to function extract_identifiable_functions_raw(
        io_equations::Dict{P, P},
        ode::ODE{P},
        known::Array{P, 1},
        with_states::Bool,
    ) where {P <: MPolyRingElem{QQFieldElem}}
    coeff_lists =
        Dict(:with_states => Array{Array{P, 1}, 1}(), :no_states => Array{Array{P, 1}, 1}())
    varnames = [var_to_str(p) for p in ode.parameters]
    if with_states
        varnames = vcat(map(var_to_str, ode.x_vars), varnames)
    end
    bring, _ = Nemo.polynomial_ring(base_ring(ode.poly_ring), varnames)

    if with_states
        @debug "Computing Lie derivatives"
        for f in states_generators(ode, io_equations)
            num, den = unpack_fraction(parent_ring_change(f, bring))
            push!(coeff_lists[:with_states], [den, num])
        end
    end

    @debug "Extracting coefficients"
    if !isempty(ode.parameters)
        nonparameters = filter(
            v -> !(var_to_str(v) in map(var_to_str, ode.parameters)),
            gens(parent(first(values(io_equations)))),
        )
        for eq in values(io_equations)
            eq_coefs = collect(values(extract_coefficients(eq, nonparameters)))
            eq_coefs = [parent_ring_change(c, bring) for c in eq_coefs]
            sort!(eq_coefs, by = (p -> (total_degree(p), length(p))))
            push!(coeff_lists[:no_states], eq_coefs)
        end
    end

    for p in known
        if all(in.(map(var_to_str, vars(p)), [map(var_to_str, gens(bring))]))
            as_list = [one(bring), parent_ring_change(p, bring)]
            if any(in.(map(var_to_str, vars(p)), [map(var_to_str, nonparameters)]))
                push!(coeff_lists[:with_states], as_list)
            else
                push!(coeff_lists[:no_states], as_list)
            end
        else
            @debug "Known quantity $p cannot be casted and is thus dropped"
        end
    end

    return coeff_lists, bring
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

The function returns a tuple containing the following:
- a list of identifiable functions (as pairs [num, denum])
- the ring containing all these functuons (either parameters only of with states)
"""
@timeit _to function initial_identifiable_functions(
        ode::ODE{T};
        prob_threshold::Float64,
        known::Array{T, 1} = Array{T, 1}(),
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

    if isempty(ode.parameters)
        @info "No parameters, so Wronskian computation is not needed"
    else
        @info "Computing Wronskians"
        flush(_si_logger[].stream)
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
    end

    add_states = (with_states == :yes) || (with_states == :maybe && !is_prime)

    id_funcs, bring = extract_identifiable_functions_raw(
        io_equations,
        ode,
        empty(ode.parameters),
        add_states,
    )

    if add_states && !isempty(ode.parameters)
        @debug "Generators of identifiable functions involve states or the ideal of IO-projections is not prime, the parameter-only part is getting simplified"
        # NOTE: switching to a ring without states for a moment
        param_ring, _ = polynomial_ring(
            base_ring(bring),
            map(string, ode.parameters),
            internal_ordering = Nemo.internal_ordering(bring),
        )
        id_funcs_no_states_param = map(
            polys -> map(poly -> parent_ring_change(poly, param_ring), polys),
            id_funcs[:no_states],
        )
        _runtime_logger[:check_time] =
            @elapsed no_states_simplified = simplified_generating_set(
            RationalFunctionField(id_funcs_no_states_param),
            prob_threshold = prob_threshold,
            seed = 42,
            simplify = :standard,
            rational_interpolator = rational_interpolator,
            cmp = cmp,
        )
        dennums_simplified = fractions_to_dennums(no_states_simplified)
        # switch back the ring
        id_funcs[:no_states] = map(
            polys -> map(poly -> parent_ring_change(poly, bring), polys),
            dennums_simplified,
        )
    end

    if !add_states
        return id_funcs[:no_states], bring
    end
    return vcat(id_funcs[:with_states], id_funcs[:no_states]), bring
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
    states_needed = false
    for f in funcs_to_check
        num, den = unpack_fraction(f)
        if !all(v -> v in ode.parameters, union(vars(num), vars(den)))
            @info "Functions to check involve states"
            states_needed = true
            break
        end
    end
    if !states_needed && isempty(ode.parameters)
        return [true for _ in funcs_to_check]
    end

    half_p = 0.5 + prob_threshold / 2
    identifiable_functions_raw, bring = initial_identifiable_functions(
        ode,
        known = known,
        prob_threshold = half_p,
        var_change_policy = var_change_policy,
        with_states = states_needed ? :yes : :maybe,
    )

    funcs_to_check = Vector{Generic.FracFieldElem{P}}(
        map(f -> parent_ring_change(f, bring) // one(bring), funcs_to_check),
    )

    _runtime_logger[:check_time] =
        get(_runtime_logger, :check_time, 0.0) + @elapsed result = field_contains(
        RationalFunctionField(identifiable_functions_raw),
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
