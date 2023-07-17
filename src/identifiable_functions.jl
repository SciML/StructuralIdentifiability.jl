
function find_identifiable_functions end

"""
    IdealMQS

Relations between elements of Q(x) encoded with an ideal in Q(x)[y]:

    <num_i(x) den_i(y) - num_i(y) den_i(x), Q(y) t - 1>

`IdealMQS` provides all of the functions that `AbstractBlackboxIdeal` requires.
See `?ParamPunPam.AbstractBlackboxIdeal` for details.

"""
mutable struct IdealMQS{PolyQQ} <: AbstractBlackboxIdeal
    funcs_den_nums::Vector{Vector{PolyQQ}}
    nums_qq::Vector{PolyQQ}
    dens_qq::Vector{PolyQQ}
    sat_qq::PolyQQ
    nums_gf::Any
    dens_gf::Any
    sat_gf::Any
    parent_ring_param::Any

    """
        IdealMQS(funcs_den_nums)
    
    ## Example

    ```jldoctest

    ```
    """
    function IdealMQS(
        funcs_den_nums::Vector{Vector{PolyQQ}};
        sat_varname = "t",
    ) where {PolyQQ}
        @assert !isempty(funcs_den_nums)
        R = parent(first(first(funcs_den_nums)))
        @info "Recording the ideal generators in $R"
        K, n, ord = base_ring(R), nvars(R), Nemo.ordering(R)
        Q = reduce(lcm, map(first, funcs_den_nums))
        @debug "Rational functions common denominator" Q
        isone(Q) && (@warn "Common denominator of the field generators is one" Q)
        existing_varnames = map(String, symbols(R))
        @info "Saturating variable if $sat_varname"
        # NOTE: what about F4-sat? 
        # NOTE: if this becomes a bottleneck, one of
        # the two ring conversions can be removed
        ystrs = ["y$i" for i in 1:length(existing_varnames)]
        R_y, _ = Nemo.PolynomialRing(K, ystrs, ordering = ord)
        funcs_den_nums = map(
            dennums ->
                map(f -> parent_ring_change(f, R_y, matching = :byindex), dennums),
            funcs_den_nums,
        )
        Q = parent_ring_change(Q, R_y, matching = :byindex)
        @assert !(sat_varname in ystrs)
        varnames = pushfirst!(ystrs, sat_varname)
        R_sat, v_sat = Nemo.PolynomialRing(K, varnames, ordering = ord)
        t_sat = first(v_sat)
        Q_sat = parent_ring_change(Q, R_sat)
        sat_qq = Q_sat * t_sat - 1
        nums_qq = empty(funcs_den_nums[1])
        dens_qq = empty(nums_qq)
        for i in 1:length(funcs_den_nums)
            dennums = funcs_den_nums[i]
            # NOTE(Alex): double-check if the generators are generated correctly 
            # @assert length(dennums) > 1 "Strange field generators: $dennums"
            den = dennums[1]
            den = parent_ring_change(den, R_sat)
            # NOTE: remove duplicates in numerators
            for j in 2:length(dennums)
                num = dennums[j]
                num = parent_ring_change(num, R_sat)
                if iszero(num)
                    continue
                end
                push!(nums_qq, num)
                push!(dens_qq, den)
            end
        end
        parent_ring_param, _ = PolynomialRing(R, varnames, ordering = ord)
        @info "Generated MQS ideal in $R_sat with $(length(nums_qq) + 1) elements"
        new{elem_type(R_sat)}(
            funcs_den_nums,
            nums_qq,
            dens_qq,
            sat_qq,
            nothing,
            nothing,
            nothing,
            parent_ring_param,
        )
    end
end

Base.length(ideal::IdealMQS) = length(ideal.nums_qq) + 1
AbstractAlgebra.base_ring(ideal::IdealMQS) = base_ring(ideal.nums_qq[1])
ParamPunPam.base_ring_mod_p(ideal::IdealMQS) = base_ring(ideal.nums_gf[1])
AbstractAlgebra.parent(ideal::IdealMQS) = ideal.parent_ring_param
ParamPunPam.parent_params(ideal::IdealMQS) = base_ring(ideal.parent_ring_param)
ParamPunPam.parent_mod_p(ideal::IdealMQS) = parent(ideal.nums_gf[1])

function ideal_generators_raw(mqs::IdealMQS)
    return field_to_ideal(mqs.funcs_den_nums)
end

function ParamPunPam.reduce_mod_p!(mqs::IdealMQS, ff)
    @info "Reducing MQS ideal modulo $(ff).."
    # TODO: check that the reduction is lucky
    nums_qq, dens_qq, sat_qq = mqs.nums_qq, mqs.dens_qq, mqs.sat_qq
    nums_gf, dens_gf = map(
        polys -> map(poly -> map_coefficients(c -> ff(c), poly), polys),
        (nums_qq, dens_qq),
    )
    sat_gf = map_coefficients(c -> ff(c), sat_qq)
    mqs.nums_gf = nums_gf
    mqs.dens_gf = dens_gf
    mqs.sat_gf = sat_gf
    return nothing
end

function ParamPunPam.evaluate_mod_p(mqs::IdealMQS, point)
    @debug "Evaluating MQS ideal at $point"
    nums_gf, dens_gf, sat_gf = mqs.nums_gf, mqs.dens_gf, mqs.sat_gf
    @assert !isnothing(nums_gf) && !isnothing(dens_gf) && !isnothing(sat_gf)
    K_1 = ParamPunPam.base_ring_mod_p(mqs)
    K_2 = parent(first(point))
    @assert K_1 == K_2
    @assert length(point) == nvars(ParamPunPam.parent_params(mqs))
    # +1 actual variable because of the saturation!
    @assert length(point) + 1 == nvars(parent(nums_gf[1]))
    # TODO: Assuming the saturating variable is the first one
    point_sat = vcat(one(K_1), point)
    nums_gf_spec = map(num -> evaluate(num, point_sat), nums_gf)
    # TODO: a single denominator can be evaluated once per batch of numerators
    dens_gf_spec = map(den -> evaluate(den, point_sat), dens_gf)
    polys = Vector{typeof(sat_gf)}(undef, length(nums_gf_spec) + 1)
    for i in 1:length(nums_gf_spec)
        num, num_spec = nums_gf[i], nums_gf_spec[i]
        den, den_spec = dens_gf[i], dens_gf_spec[i]
        polys[i] = num * den_spec - den * num_spec
    end
    polys[end] = sat_gf
    return polys
end

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

"""
    dennums_to_fractions(dennums)
    
Returns the field generators represented by fractions.

Input: an array of arrays of polynomials, as in 
`[[f1, f2, f3, ...], [g1, g2, g3, ...], ...]`

Output: an array of fractions `[f2/f1, f3/f1, ..., g2/g1, g3/g1, ...]`
"""
function dennums_to_fractions(dennums::Vector{Vector{T}}) where {T}
    fractions = Vector{AbstractAlgebra.Generic.Frac{T}}()
    for dni in dennums
        den, nums = dni[1], dni[2:end]
        append!(fractions, map(c -> c // den, nums))
    end
    fractions
end

"""
    fractions_to_dennums(fractions)
    
Returns the field generators represented by denominators and numerators.

Input: an array of fractions, as in `[f2/f1, f3/f1, ..., g2/g1, g3/g1, ...]`

Output: an array of arrays of polynomials, as in 
`[[f1, f2, f3, ...], [g1, g2, g3, ...], ...]`
"""
function fractions_to_dennums(fractions)
    # NOTE: Maybe collapse similar denominators
    map(f -> [denominator(f), numerator(f)], fractions)
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

function field_generators(ode::ODE{T}, p::Float64 = 0.99) where {T}
    @info "Computing IO-equations"
    runtime = @elapsed io_equations = find_ioequations(ode)
    @info "IO-equations computed in $runtime seconds"
    _runtime_logger[:id_io_time] = runtime

    @info "Assessing global identifiability"
    runtime = @elapsed global_result = check_identifiability(
        io_equations,
        ode,
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
    funcs_den_nums
end

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
    Rt, vt = Nemo.PolynomialRing(K, varnames, ordering = ord)
    if append_front
        xs, t = vt[2:end], first(vt)
    else
        xs, t = vt[1:(end - 1)], last(vt)
    end
    It = map(f -> parent_ring_change(f, Rt), I)
    Qt = parent_ring_change(Q, Rt)
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
    @info "Producing the ideal generators in $R"
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
    Fy = map(f -> parent_ring_change(f, Ry, matching = :byindex), Fx)
    Qy = parent_ring_change(Qx, Ry, matching = :byindex)
    # NOTE(Alex): I think we want to clean up Fx before creating the
    # ideal. For example, it is possible that some of Fx are duplicates. Would
    # make sense to do so in ParamPunPam.jl, but here we have an advantage that
    # the polynomials are not constructed yet
    #
    # one backward gaussian elimination pass for Fx?
    #
    I = empty(ys)
    for (Fyi, Fxi) in zip(Fy, Fx)
        Gx = gcd(Qx, Fxi)
        Qxg = divexact(Qx, Gx)
        Fxig = divexact(Fxi, Gx)
        F = Fyi * Qxg - Fxig * Qy
        push!(I, F)
    end
    I, t = saturate(I, Qy)
    # TODO: remove this line
    I_rat = map(f -> map_coefficients(c -> c // one(R), f), I)
    return I_rat
end

"""
    simplify_identifiable_functions(X)

TODO: describe the format of the output.
"""
function simplify_identifiable_functions(
    funcs_den_nums::Vector{Vector{T}},
    p = 0.99,
    seed = 42,
) where {T}
    # TODO: use seed!
    # Generate the ideal
    @info "Simplifying identifiable functions"
    runtime = @elapsed mqs = IdealMQS(funcs_den_nums)
    @info "Differential ideal created in $runtime seconds"
    _runtime_logger[:id_ideal_time] = runtime
    # Find a GB
    # NOTE(Alex): I think correctness checks, tweaks with the total degree of
    # interpolation, and other GB implementation details belong in
    # ParamPunPam.jl
    # NOTE(Alex): The algorithms highly depend on randomness, we really want to
    # fix the random seed here. Maybe provide a keyword argument in
    # find_identifiable_functions
    @info "Computing a Groebner basis"
    _runtime_logger[:id_groebner_time] = 0.0
    two_sided_inclusion = false
    gb = nothing
    current_degrees = (2, 2)
    while !two_sided_inclusion
        @info "Computing GB with parameters up to degrees $(current_degrees)"
        runtime = @elapsed gb = ParamPunPam.paramgb(mqs, up_to_degree = current_degrees)
        _runtime_logger[:id_groebner_time] += runtime
        @info "Groebner basis computed in $runtime seconds"

        id_coeffs = map(collect ∘ coefficients, gb)
        id_coeffs_set = Set{AbstractAlgebra.Generic.Frac{T}}()
        for id_coeff in id_coeffs
            union!(id_coeffs_set, id_coeff)
        end
        id_funcs = collect(id_coeffs_set)
        @info "Identifiable functions up to degrees $(current_degrees) are" id_funcs
        @info "Checking two-sided inclusion with probability $p"
        time_start = time_ns()
        id_funcs_den_nums = fractions_to_dennums(id_funcs)
        original_id_funcs = dennums_to_fractions(funcs_den_nums)
        # Check inclusion: <simplified generators> in <original generators> 
        inclusion = check_field_membership(id_funcs_den_nums, original_id_funcs, p)
        two_sided_inclusion = two_sided_inclusion || all(inclusion)
        # Check inclusion: <original generators> in <simplified generators>
        inclusion = check_field_membership(funcs_den_nums, id_funcs, p)
        _runtime_logger[:id_inclusion_check] = (time_ns() - time_start) / 1e9
        @info "Inclusion checked in $(_runtime_logger[:id_inclusion_check]) seconds"
        two_sided_inclusion = two_sided_inclusion && all(inclusion)
        current_degrees = current_degrees .* 2
    end

    # Simplify.
    # NOTE(Alex): We can check redundancy of each of the new generators, but
    # let's not do that for now.
    id_coeffs = map(collect ∘ coefficients, gb)
    id_coeffs_set = Set{AbstractAlgebra.Generic.Frac{T}}()
    for id_coeff in id_coeffs
        union!(id_coeffs_set, id_coeff)
    end
    @info "The coefficients of the Groebner basis are presented by $(length(id_coeffs_set)) rational functions"
    time_start = time_ns()
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
    _runtime_logger[:id_filter_time] = (time_ns() - time_start) / 1e9
    @info "Functions filtered and sorted in $(_runtime_logger[:id_filter_time]) seconds"
    @info "$(length(id_funcs)) functions after simplification"
    # Convert back into the [denominator, numerators...] format
    id_funcs_den_nums = fractions_to_dennums(id_funcs)
    isempty(id_funcs) && return id_funcs_den_nums
    id_funcs_den_nums
end

# TODO: add a meaningful simplification step for field generators
"""
    find_identifiable_functions(ode::ODE, p = 0.99)

Finds all functions of parameters that are identifiable for the given ODE
system.

Input:
- `ode` - `ODE`-system.
Output:
- a set of generators of the field of identifiable functions. Any function in
  the returned set (or a combination thereof) is globally identifiable.

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
    ode::ODE{T},
    p::Float64 = 0.99;
    simplify = true,
    seed = 42,
) where {T <: MPolyElem{fmpq}}
    runtime_start = time_ns()

    # TODO: Can we still provide parameter `p=0.99` here?
    runtime = @elapsed funcs_den_nums = field_generators(ode, p)
    # TODO: Check if all parameters are globally id.
    _runtime_logger[:id_field_generators] = runtime

    if simplify
        funcs_den_nums = simplify_identifiable_functions(funcs_den_nums, p, seed)
    end
    id_funcs = dennums_to_fractions(funcs_den_nums)
    _runtime_logger[:id_total] = (time_ns() - runtime_start) / 1e9
    @info "The search for identifiable functions concluded in $(_runtime_logger[:id_total]) seconds"
    return id_funcs
end

"""
    find_identifiable_functions(ode::ModelingToolkit.ODESystem; measured_quantities=Array{ModelingToolkit.Equation}[])

Finds all functions of parameters that are identifiable for the given ODE
system.
    
Input:
- `ode` - the ModelingToolkit.ODESystem object that defines the model
- `measured_quantities` - the output functions of the model

"""
function find_identifiable_functions(
    ode::ModelingToolkit.ODESystem;
    measured_quantities = Array{ModelingToolkit.Equation}[],
    p = 0.99,
    simplify = true,
    seed = 42,
)
    if length(measured_quantities) == 0
        if any(ModelingToolkit.isoutput(eq.lhs) for eq in ModelingToolkit.equations(ode))
            @info "Measured quantities are not provided, trying to find the outputs in input ODE."
            measured_quantities = filter(
                eq -> (ModelingToolkit.isoutput(eq.lhs)),
                ModelingToolkit.equations(ode),
            )
        else
            throw(
                error(
                    "Measured quantities (output functions) were not provided and no outputs were found.",
                ),
            )
        end
    end
    params = ModelingToolkit.parameters(ode)
    ode, conversion = preprocess_ode(ode, measured_quantities)
    out_funcs = Vector{Num}()
    params_ = [eval_at_nemo(each, conversion) for each in params]
    result = find_identifiable_functions(ode, p, simplify = simplify, seed = seed)
    nemo2mtk = Dict(params_ .=> map(Num, params))
    out_funcs = [eval_at_dict(func, nemo2mtk) for func in result]
    return out_funcs
end
