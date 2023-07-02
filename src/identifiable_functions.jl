
"""
    extract_identifiable_functions(io_equations, parameters)
For the io_equation and the list of all parameter variables, returns a set of generators of a field of all functions of parameters
Note: an experimental functionality at the moment, may fail and may be inefficient
"""
function extract_identifiable_functions(
    io_equations::Array{P, 1},
    parameters::Array{P, 1},
    known_functions::Array{P, 1},
) where {P <: MPolyElem{fmpq}}
    @debug "Extracting coefficients"
    flush(stdout)
    nonparameters = filter(
        v -> !(var_to_str(v) in map(var_to_str, parameters)),
        gens(parent(io_equations[1])),
    )
    coeff_lists = Array{Array{P, 1}, 1}()
    for eq in io_equations
        push!(coeff_lists, collect(values(extract_coefficients(eq, nonparameters))))
    end
    R = parent(first(first(coeff_lists)))
    for f in known_functions
        push!(coeff_lists, [one(R), parent_ring_change(f, R)])
    end
    for p in coeff_lists
        @debug sort(map(total_degree, p))
    end

    @debug "Resulting Coefficient List: $coeff_lists"

    return coeff_lists
end

#------------------------------------------------------------------------------

"""
    extract_identifiable_functions_raw(io_equations, parameters)
For the io_equation and the list of all parameter variables, returns a set of *raw* *generators of a field of all functions of parameters
"""
function extract_identifiable_functions_raw(
    io_equations::Array{P, 1},
    parameters::Array{P, 1},
) where {P <: MPolyElem{fmpq}}
    @debug "Extracting coefficients"
    flush(stdout)
    nonparameters = filter(
        v -> !(var_to_str(v) in map(var_to_str, parameters)),
        gens(parent(io_equations[1])),
    )
    result = []
    for eq in io_equations
        coeffs = sort(
            collect(values(extract_coefficients(eq, nonparameters))),
            by = total_degree,
        )
        append!(result, [c // first(coeffs) for c in coeffs[2:end]])
    end

    return result
end

#------------------------------------------------------------------------------

"""
    saturate(I, Q)

"""
function saturate(I, Q; varname = "t", append_front = true)
    @info "Saturating the ideal, saturating variable is $varname"
    R = parent(first(I))
    K, ord = base_ring(R), ordering(R)
    existing_varnames = map(String, symbols(R))
    @assert !(varname in existing_varnames)
    varnames =
        append_front ? pushfirst!(existing_varnames, varname) :
        push!(existing_varnames, varname)
    _, vt = Nemo.PolynomialRing(K, varnames, ordering = ord)
    if append_front
        xs, t = vt[2:end], first(vt)
    else
        xs, t = vt[1:(end - 1)], last(vt)
    end
    It = map(f -> evaluate(f, xs), I)
    Qt = evaluate(Q, xs)
    sat = 1 - Qt * t
    push!(It, sat)
    It, t
end

"""
    field_to_ideal(X)

"""
function field_to_ideal(
    funcs_den_nums::Vector{Vector{T}};
    top_level_var = "y",
    top_level_ord = :degrevlex,
) where {T}
    @assert !isempty(funcs_den_nums)
    R = parent(first(first(funcs_den_nums)))
    K, n = base_ring(R), nvars(R)
    Q = reduce(lcm, map(first, funcs_den_nums))
    @debug "Rational functions common denominator" Q
    Fx, Qx = Vector{T}(), Q
    for den_nums in funcs_den_nums
        den, nums = den_nums[1], den_nums[2:end]
        sep = divexact(Q, den)
        append!(Fx, nums .* sep)
    end
    ystrs = ["$top_level_var$i" for i in 1:n]
    Ry, ys = Nemo.PolynomialRing(R, ystrs, ordering = top_level_ord)
    Fy = map(f -> change_base_ring(R, f, parent = Ry), Fx)
    Qy = change_base_ring(R, Q, parent = Ry)
    # TODO: ensure that the polys are primitive
    I = [Fyi * Qx - Fxi * Qy for (Fyi, Fxi) in zip(Fy, Fx)]
    I, t = saturate(I, Qy)
    # TODO: remove this line
    I_rat = map(f -> map_coefficients(c -> c // one(R), f), I)
    return I_rat
end

is_ratfunc_const(f) = is_constant(numerator(f)) && is_constant(denominator(f))
is_ratfunc_normalized(f) =
    (leading_coefficient(denominator(f)) > 0) && isone(gcd(numerator(f), denominator(f)))

"""
    ratfunc_cmp(f, g)

Returns `f < g`.

`f < g` iff the total degree of the numerator of `f` is less than that of `g`.
Breaks ties by comparing num(f) < num(g) and den(f) < den(g) with the current
monomial ordering.
"""
function ratfunc_cmp(f, g)
    numf, denf = unpack_fraction(f)
    numg, deng = unpack_fraction(g)
    total_degree(numf) < total_degree(numg) && return true
    total_degree(numf) > total_degree(numg) && return false
    leading_monomial(numf) < leading_monomial(numg) && return true
    leading_monomial(numf) > leading_monomial(numg) && return false
    leading_monomial(denf) < leading_monomial(deng) && return true
    return false
end

"""
    simplify_identifiable_functions(X)

TODO: describe the format of the output.
"""
function simplify_identifiable_functions(funcs_den_nums::Vector{Vector{T}}) where {T}
    # Generate the ideal
    @info "Simplifying identifiable functions"
    runtime = @elapsed dideal = field_to_ideal(funcs_den_nums)
    @info "Differential ideal computed in $runtime seconds"
    _runtime_logger[:id_ideal_time] = runtime
    # Find a GB
    # NOTE(Alex): I think correctness checks, tweaks with the total degree of
    # interpolation, and other GB implementation details belong in
    # ParamPunPam.jl
    # NOTE(Alex): The algorithms highly depend on randomness, we really want to
    # fix the random seed here. Maybe provide a keyword argument in
    # find_identifiable_functions
    @info "Computing a Groebner basis"
    @debug "The polynomial ring is $(parent(first(dideal)))"
    runtime = @elapsed gb = ParamPunPam.paramgb(dideal)
    @info "Groebner basis computed in $runtime seconds"
    _runtime_logger[:id_groebner_time] = runtime
    # Simplify.
    # NOTE(Alex): We can check redundancy of each of the new generators, but
    # let's not do that for now.
    id_coeffs = map(collect ∘ coefficients, gb)
    id_coeffs_set = Set{AbstractAlgebra.Generic.Frac{T}}()
    for id_coeff in id_coeffs
        union!(id_coeffs_set, id_coeff)
    end
    @info "The coefficients of the Groebner basis are presented by $(length(id_coeffs_set)) rational functions"
    id_funcs = collect(id_coeffs_set)
    filter!(!is_ratfunc_const, id_funcs)
    @assert all(is_ratfunc_normalized, id_funcs)
    for i in 1:length(id_funcs)
        func = id_funcs[i]
        num, den = unpack_fraction(func)
        if is_constant(num)
            func = den // num
        end
        # asserting they are comparable!
        lcn = leading_coefficient(num)
        if lcn < 0
            func = func * lcn
        end
        id_funcs[i] = func
    end
    sort!(id_funcs, lt = ratfunc_cmp)
    @info "$(length(id_funcs)) functions after simplification"
    # Convert back into the [denominator, numerators...] format
    id_funcs_den_nums = map(f -> [denominator(f), numerator(f)], id_funcs)
    # Check inclusion: <original generators> in <simplified generators>
    # TODO:
    id_funcs_den_nums
end

#------------------------------------------------------------------------------

"""
    find_identifiable_functions(ode::ODE)

Finds all functions of parameters that are identifiable for the given ODE
system.

Input:
- `ode` - `ODE`-system.

Output:
- a set of generators of the field of identifiable functions. Any function in
the returned set (or a combination thereof) is globally identifiable.
"""
function find_identifiable_functions(
    ode::ODE{T},
    p::Float64 = 0.99;
    simplify = true,
) where {T <: MPolyElem{fmpq}}
    # TODO: Can we still provide parameter `p=0.99` here?
    runtime_start = time_ns()
    @info "Computing IO-equations"
    runtime = @elapsed io_equations = find_ioequations(ode)
    @info "IO-equations computed in $runtime seconds"
    _runtime_logger[:id_io_time] = runtime

    @info "Assessing global identifiability"
    runtime = @elapsed global_result = check_identifiability(
        collect(values(io_equations)),
        ode.parameters,
        empty(ode.parameters),
        p,
    )
    @info "Global identifiability assessed in $runtime seconds"
    _runtime_logger[:id_global_time] = runtime

    known_params = Array{fmpq_mpoly, 1}()
    for (glob, p) in zip(global_result, ode.parameters)
        if glob
            push!(known_params, p)
        end
    end
    runtime = @elapsed funcs_den_nums = extract_identifiable_functions(
        collect(values(io_equations)),
        ode.parameters,
        known_params,
    )
    _runtime_logger[:id_extract_funcs_time] = runtime
    # NOTE(Alex): perhaps we would want to get rid of the back and forth
    # conversion of formats at some point
    if simplify
        funcs_den_nums = simplify_identifiable_functions(funcs_den_nums)
    end
    id_funcs = Vector{AbstractAlgebra.Generic.Frac{T}}()
    for den_nums in funcs_den_nums
        den, nums = den_nums[1], den_nums[2:end]
        append!(id_funcs, map(c -> c // den, nums))
    end
    _runtime_logger[:id_total] = time_ns() - runtime_start
    @info "The search for identifiable functions concluded in $(_runtime_logger[:id_total]) seconds"
    return id_funcs
end
