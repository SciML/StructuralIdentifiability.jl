
function extract_identifiable_functions_raw(
    io_equations::Dict{P, P},
    ode::ODE{P},
    known::Array{P, 1},
    with_states::Bool,
) where {P <: MPolyElem{fmpq}}
    coeff_lists = Array{Array{P, 1}, 1}()
    bring = nothing
    if with_states
        @debug "Computing Lie derivatives"
        for f in states_generators(ode, io_equations)
            num, den = unpack_fraction(f)
            push!(coeff_lists, [den, num])
        end
        bring = parent(first(first(coeff_lists)))
    end

    @debug "Extracting coefficients"
    flush(stdout)
    nonparameters = filter(
        v -> !(var_to_str(v) in map(var_to_str, ode.parameters)),
        gens(parent(first(values(io_equations)))),
    )
    for eq in values(io_equations)
        eq_coefs = collect(values(extract_coefficients(eq, nonparameters)))
        if !isnothing(bring)
            eq_coefs = [parent_ring_change(c, bring) for c in eq_coefs]
        end
        push!(coeff_lists, eq_coefs)
    end
    if isnothing(bring)
        bring = parent(first(first(coeff_lists)))
    end
    for p in known
        if all(in.(map(var_to_str, vars(p)), [map(var_to_str, gens(bring))]))
            push!(coeff_lists, [one(bring), parent_ring_change(p, bring)])
        else
            @debug "Known quantity $p cannot be casted and is thus dropped"
        end
    end

    # NOTE: Returned entities live in a new ring, different from the one
    # attached to the input ODE. 
    # The new ring includes only parameter variables (and optionally states).
    new_vars = ode.parameters
    if with_states
        new_vars = vcat(new_vars, ode.x_vars)
    end
    new_ring, _ = PolynomialRing(Nemo.QQ, map(Symbol, new_vars))
    coeff_lists =
        map(coeffs -> map(c -> parent_ring_change(c, new_ring), coeffs), coeff_lists)

    return coeff_lists
end

# ------------------------------------------------------------------------------

function check_identifiability(
    io_equations::Dict{P, P},
    ode::ODE{P},
    known::Array{P, 1},
    funcs_to_check::Array{<:Any, 1},
    p::Float64 = 0.99,
) where {P <: MPolyElem{fmpq}}
    states_needed = false
    for f in funcs_to_check
        num, den = unpack_fraction(f)
        if !all(v -> v in ode.parameters, union(vars(num), vars(den)))
            @info "Functions to check involve states"
            states_needed = true
            break
        end
    end

    identifiable_functions_raw =
        extract_identifiable_functions_raw(io_equations, ode, known, states_needed)
    bring = parent(first(first(identifiable_functions_raw)))

    funcs_to_check = Vector{Generic.Frac{P}}(
        map(f -> parent_ring_change(f, bring) // one(bring), funcs_to_check),
    )

    return field_contains(
        RationalFunctionField(identifiable_functions_raw),
        funcs_to_check,
        p,
    )
end

function check_identifiability(
    io_equations::Dict{P, P},
    ode::ODE{P},
    known::Array{P, 1},
    p::Float64 = 0.99,
) where {P <: MPolyElem{fmpq}}
    return check_identifiability(io_equations, ode, known, ode.parameters, p)
end

#------------------------------------------------------------------------------
"""
    assess_global_identifiability(ode::ODE{P}, p::Float64=0.99; var_change=:default) where P <: MPolyElem{fmpq}

Input:
- `ode` - the ODE model
- `known` - a list of functions in states which are assumed to be known and generic
- `p` - probability of correctness
- `var_change` - a policy for variable change (`:default`, `:yes`, `:no`), affects only the runtime

Output:
- a dictionary mapping each parameter to a boolean.

Checks global identifiability for parameters of the model provided in `ode`. Call this function to check global identifiability of all parameters automatically.
"""
function assess_global_identifiability(
    ode::ODE{P},
    known::Array{P, 1} = Array{P, 1}(),
    p::Float64 = 0.99;
    var_change = :default,
) where {P <: MPolyElem{fmpq}}
    result_list = assess_global_identifiability(
        ode,
        ode.parameters,
        known,
        p;
        var_change = var_change,
    )

    return Dict(param => val for (param, val) in zip(ode.parameters, result_list))
end

#------------------------------------------------------------------------------

"""
    assess_global_identifiability(ode, [funcs_to_check, p=0.99, var_change=:default])

Input:
- `ode` - the ODE model
- `funcs_to_check` - rational functions in parameters
- `known` - function in parameters that are assumed to be known and generic
- `p` - probability of correctness
- `var_change` - a policy for variable change (`:default`, `:yes`, `:no`),
                affects only the runtime

Output:
- array of length `length(funcs_to_check)` with true/false values for global identifiability
        or dictionary `param => Bool` if `funcs_to_check` are not given

Checks global identifiability of functions of parameters specified in `funcs_to_check`.
"""
function assess_global_identifiability(
    ode::ODE{P},
    funcs_to_check::Array{<:Any, 1},
    known::Array{P, 1} = Array{P, 1}(),
    p::Float64 = 0.99;
    var_change = :default,
) where {P <: MPolyElem{fmpq}}
    submodels = find_submodels(ode)
    if length(submodels) > 0
        @info "Note: the input model has nontrivial submodels. If the computation for the full model will be too heavy, you may want to try to first analyze one of the submodels. They can be produced using function `find_submodels`"
    end

    @info "Computing IO-equations"
    ioeq_time =
        @elapsed io_equations = find_ioequations(ode; var_change_policy = var_change)
    @debug "Sizes: $(map(length, values(io_equations)))"
    @info "Computed in $ioeq_time seconds" :ioeq_time ioeq_time
    _runtime_logger[:ioeq_time] = ioeq_time

    @info "Computing Wronskians"
    wrnsk_time = @elapsed wrnsk = wronskian(io_equations, ode)
    @info "Computed in $wrnsk_time seconds" :wrnsk_time wrnsk_time
    _runtime_logger[:wrnsk_time] = wrnsk_time

    dims = map(ncols, wrnsk)
    @info "Dimensions of the Wronskians $dims"

    rank_times = @elapsed wranks = map(rank, wrnsk)
    @debug "Dimensions of the Wronskians $dims"
    @debug "Ranks of the Wronskians $wranks"
    @info "Ranks of the Wronskians computed in $rank_times seconds" :rank_time rank_times
    _runtime_logger[:rank_time] = rank_times

    if any([dim != rk + 1 for (dim, rk) in zip(dims, wranks)])
        @warn "One of the Wronskians has corank greater than one, so the results of the algorithm will be valid only for multiexperiment identifiability. If you still  would like to assess single-experiment identifiability, we recommend using SIAN (https://github.com/alexeyovchinnikov/SIAN-Julia)"
    end

    @info "Assessing global identifiability using the coefficients of the io-equations"
    check_time =
        @elapsed result = check_identifiability(io_equations, ode, known, funcs_to_check, p)
    @info "Computed in $check_time seconds" :check_time check_time
    _runtime_logger[:check_time] = check_time

    return result
end
