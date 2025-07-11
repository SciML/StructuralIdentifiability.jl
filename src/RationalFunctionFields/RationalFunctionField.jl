
"""
    RationalFunctionField

A subfield of the field of rational functions over the rationals.

## Example

```jldoctest
using Nemo
using StructuralIdentifiability: RationalFunctionField

R, (x, y, z) = QQ["x", "y", "z"]

# Constructs a subfield generated by x / y, y / z
rff = RationalFunctionField([x // y, y // z])

# Constructs a subfield generated by y / x, 1 / x, z / y
rff = RationalFunctionField([[x, y, R(1)], [y, z]])
```
"""
mutable struct RationalFunctionField{T}
    dennums::Vector{Vector{T}}
    mqs::IdealMQS{T}
    mqs_membership::IdealMQS{T}

    # cached transcendence-related information
    trbasis_probability::Float64
    trbasis::Vector{Generic.FracFieldElem{T}}
    # trancendence basis of the ambient rational function field over the given one
    trbasis_over::Vector{T}

    function RationalFunctionField(polys::Vector{T}) where {T}
        RationalFunctionField(polys .// one(parent(first(polys))))
    end
    function RationalFunctionField(fractions::Vector{Generic.FracFieldElem{T}}) where {T}
        RationalFunctionField(fractions_to_dennums(fractions))
    end
    function RationalFunctionField(dennums::Vector{Vector{T}}) where {T}
        @assert !isempty(dennums)
        F = new{T}(
            dennums,
            IdealMQS(dennums),
            IdealMQS(dennums),
            0,
            Vector{Generic.FracFieldElem{T}}(),
            Vector{T}(),
        )
        update_trbasis_info!(F, 0.9999)
        return F
    end
end

# ------------------------------------------------------------------------------

function trivial(F::RationalFunctionField)
    return isempty(F.dennums)
end

function poly_ring(F::RationalFunctionField)
    return parent(first(first(F.dennums)))
end

function generators(F::RationalFunctionField)
    return dennums_to_fractions(F.dennums)
end

function Base.zero(F::RationalFunctionField)
    return zero(poly_ring(F)) // one(poly_ring(F))
end

function Base.one(F::RationalFunctionField)
    return one(poly_ring(F)) // one(poly_ring(F))
end

# ------------------------------------------------------------------------------

function update_trbasis_info!(F::RationalFunctionField, p::Float64)
    F.trbasis_probability = p
    fgens = generators(F)
    base_vars = gens(poly_ring(F))
    if isempty(base_vars)
        return
    end
    if isempty(fgens)
        return
    end
    maxdeg = maximum(map(total_degree_frac, fgens), init = 1) - 1
    # degree of the polynomial whose nonvanishing will be needed for correct result
    D = max(10, Int(ceil(maxdeg * length(base_vars) / (1 - p))))
    eval_point = [Nemo.QQ(rand(1:D)) for x in base_vars]

    J = jacobian(fgens, eval_point)
    pivots, _ = select_pivots(Nemo.rref(J)[2])
    _, nonpivots = select_pivots(Nemo.rref(transpose(J))[2])

    old_trbasis = F.trbasis
    F.trbasis = [fgens[i] for i in pivots]
    F.trbasis_over = [base_vars[i] for i in nonpivots]
    @assert length(F.trbasis) + length(F.trbasis_over) == length(base_vars)

    if old_trbasis != F.trbasis
        F.mqs_membership =
            IdealMQS(vcat(F.dennums, [[x, one(poly_ring(F))] for x in F.trbasis_over]))
    end
end

# ------------------------------------------------------------------------------

function _check_algebraicity(trbasis, ratfuncs, sampling_range)
    if isempty(ratfuncs)
        return Bool[]
    end
    if isempty(trbasis)
        return map(f -> total_degree_frac(f) == 0, ratfuncs)
    end
    polyring = parent(numerator(first(trbasis)))
    field = base_ring(polyring)

    while true
        eval_point = [field(rand(1:sampling_range)) for _ in gens(polyring)]

        J = jacobian(vcat(trbasis, zero(first(trbasis))), eval_point)
        rank = LinearAlgebra.rank(J)
        if rank < length(trbasis)
            continue
        end

        result = Bool[]
        for f in ratfuncs
            f = parent_ring_change(f, polyring)
            for (j, x) in enumerate(gens(polyring))
                J[j, end] = evaluate(derivative(f, x), eval_point)
            end
            push!(result, LinearAlgebra.rank(J) == rank)
        end
        return result
    end
end

# ------------------------------------------------------------------------------

"""
    check_algebraicity(field, ratfuncs, p)

Checks whether given rational function `ratfuncs` are algebraic over the field `field`
The result is correct with probability at least `p`

Inputs:
- `F` - a rational function field
- `ratfuncs` - a list of lists of rational functions.
- `p` real number from (0, 1)

Output:
- a list `L[i]` of bools of length `length(rat_funcs)` such that `L[i]` is true iff
   the i-th function is algebraic over the `field`
"""
function check_algebraicity(F::RationalFunctionField, ratfuncs, p)
    if isempty(ratfuncs)
        return Bool[]
    end
    if p > F.trbasis_probability
        update_trbasis_info!(F, p)
    end
    trbasis = F.trbasis
    maxdeg = maximum(map(total_degree_frac, vcat(ratfuncs, trbasis))) - 1

    # Here the story for correctness is tricky. Consider the cases when the answer may be wrong
    # - if the element is algebraic, then the only way the function would return incorrect result
    #   is if the trbasis was incorrect
    # - if the element is transcendental, the only way to return incorrect result would be
    #   to have the determinant vanishing regardless of the correctness of the transcendence basis

    # degree of the polynomial whose nonvanishing will be needed for correct result
    D = max(
        10,
        Int(ceil(maxdeg * (length(trbasis) + 1) * (length(ratfuncs) + 1) / (1 - p))),
    )

    return _check_algebraicity(trbasis, ratfuncs, D)
end

# ------------------------------------------------------------------------------

"""
    check_algebraicity_modp(field, ratfuncs, prime)

Checks whether given rational function `ratfuncs` are algebraic over the field `field`
via randomization modulo the given `prime`

Inputs:
- `F` - a rational function field
- `ratfuncs` - a list of lists of rational functions.
- `prime` a prime number

Output:
- a list `L[i]` of bools of length `length(rat_funcs)` such that `L[i]` is true iff
  the modular test concludes that the i-th function is algebraic over the `field`
  (no mathematical guarantees)
"""
function check_algebraicity_modp(F::RationalFunctionField, ratfuncs, prime = 2^31 - 1)
    if isempty(ratfuncs)
        return Bool[]
    end
    finite_field = Nemo.Native.GF(prime)
    trbasis_modp = [_reduce_mod_p(f, prime) for f in F.trbasis]
    ratfuncs_modp = [_reduce_mod_p(f, prime) for f in ratfuncs]

    return _check_algebraicity(trbasis_modp, ratfuncs_modp, prime - 1)
end

# ------------------------------------------------------------------------------

"""
    field_contains_mod_p(field, rat_funcs, prime)

Checks whether given rational functions belong to a given field of rational
functions over integers via a reduction modulo a prime (thus, no guarantees)
Inputs:
- `field` - a rational function field
- `ratfuncs` - a list of rational functions
- `prime` - a prime number

Output:
- a list `L[i]` of bools of length `length(rat_funcs)` such that `L[i]` is true iff
   the i-th function belongs to ``field`` (indicated by the mod-p test)
"""

@timeit _to function field_contains_mod_p(
    field::RationalFunctionField{T},
    ratfuncs::Vector{Generic.FracFieldElem{T}},
    prime = 2^31 - 1,
) where {T}
    if isempty(ratfuncs)
        return Bool[]
    end

    algebraicity = check_algebraicity_modp(field, ratfuncs, prime)
    if !any(algebraicity)
        return algebraicity
    end
    ratfuncs_algebraic = ratfuncs[algebraicity]

    ff = Nemo.Native.GF(prime)
    mqs_generators = field.mqs_membership
    reduce_mod_p!(mqs_generators, ff)

    param_ring = ParamPunPam.parent_params(mqs_generators)
    point = ParamPunPam.distinct_nonzero_points(ff, nvars(param_ring))

    gens_specialized = ParamPunPam.specialize_mod_p(mqs_generators, point)
    ratfuncs_mqs_specialized =
        specialize_fracs_to_mqs(mqs_generators, ratfuncs_algebraic, point)
    @assert parent(first(gens_specialized)) == parent(first(ratfuncs_mqs_specialized))
    gb = groebner(gens_specialized)
    d_nf = normalform(gb, ratfuncs_mqs_specialized)
    result = map(iszero, d_nf)
    return merge_results(algebraicity, result)
end

function field_contains_mod_p(
    field::RationalFunctionField{T},
    ratfuncs::Vector{Vector{T}},
    prime = 2^31 - 1,
) where {T}
    return field_contains_mod_p(field, dennums_to_fractions(ratfuncs), prime)
end

function issubfield_mod_p(
    F::RationalFunctionField{T},
    E::RationalFunctionField{T},
    prime = 2^31 - 1,
) where {T}
    return all(field_contains_mod_p(E, F.dennums, prime))
end

# ------------------------------------------------------------------------------

"""
    field_contains(field, ratfuncs, prob_threshold)

Checks whether given rational function field `field` contains given rational
functions `ratfuncs`. The result is correct with probability at least `prob_threshold`

Inputs:
- `field` - a rational function field
- `ratfuncs` - a list of rational functions
- `prob_threshold` real number from (0, 1)

Output:
- a list `L[i]` of bools of length `length(rat_funcs)` such that `L[i]` is true iff
   the i-th function belongs to `field`
"""
@timeit _to function field_contains(
    field::RationalFunctionField{T},
    ratfuncs::Vector{Generic.FracFieldElem{T}},
    prob_threshold,
) where {T}
    if isempty(ratfuncs)
        return Bool[]
    end

    half_p = 1 - (1 - prob_threshold) / 2

    algebraicity = check_algebraicity(field, ratfuncs, half_p)
    ratfuncs_algebraic = ratfuncs[algebraicity]
    if isempty(ratfuncs_algebraic)
        return algebraicity
    end

    @debug "Estimating the sampling bound"

    # uses Theorem 3.3 from https://arxiv.org/pdf/2111.00991.pdf
    # the comments below use the notation from the theorem
    ratfuncs_algebraic = [
        (iszero(f) || (total_degree(numerator(f)) > total_degree(denominator(f)))) ? f :
        1 // f for f in ratfuncs_algebraic
    ]
    denoms = map(denominator, ratfuncs_algebraic)
    ring = parent(numerator(first(ratfuncs_algebraic)))
    den_lcm = lcm(field.mqs.den_lcm_orig, foldl(lcm, denoms))
    @debug "Common lcm is $den_lcm"

    # this is deg(g) + 1
    degree = total_degree(den_lcm) + 1
    # computing maximum of deg(f) for different f's to be tested
    for (i, f) in enumerate(ratfuncs_algebraic)
        extra_degree = total_degree(den_lcm) - total_degree(denominator(f))
        degree = max(degree, extra_degree + total_degree(numerator(f)))
    end
    # computing maximum of deg(f_i) for the generators of the field
    for (i, plist) in enumerate(field.dennums)
        extra_degree = total_degree(den_lcm) - total_degree(field.mqs.dens_qq[i])
        degree = max(degree, extra_degree + maximum(total_degree, plist))
    end
    @debug "\tBound for the degrees is $degree"

    total_vars = foldl(
        union,
        map(plist -> foldl(union, map(poly -> Set(vars(poly)), plist)), field.dennums),
    )
    @debug "\tThe total number of variables in $(length(total_vars))"
    sampling_bound = BigInt(
        3 *
        BigInt(degree)^(length(total_vars) + 3) *
        (length(ratfuncs_algebraic)) *
        ceil(1 / (1 - prob_threshold)),
    )

    @debug "\tSampling from $(-sampling_bound) to $(sampling_bound)"
    mqs = field.mqs_membership
    param_ring = ParamPunPam.parent_params(mqs)
    point = map(v -> Nemo.QQ(rand((-sampling_bound):sampling_bound)), gens(param_ring))
    mqs_specialized = specialize(mqs, point)
    @debug "Computing Groebner basis ($(length(mqs_specialized)) equations)"
    mqs_ratfuncs = specialize_fracs_to_mqs(mqs, ratfuncs_algebraic, point)
    @assert parent(first(mqs_specialized)) == parent(first(mqs_ratfuncs))
    @debug "Starting the groebner basis computation"
    gb = groebner(mqs_specialized)
    result_global = map(iszero, normalform(gb, mqs_ratfuncs))
    return merge_results(algebraicity, result_global)
end

function field_contains(
    field::RationalFunctionField{T},
    ratfuncs::Vector{Vector{T}},
    prob_threshold,
) where {T}
    return field_contains(field, dennums_to_fractions(ratfuncs), prob_threshold)
end

function field_contains(
    field::RationalFunctionField{T},
    polys::Vector{T},
    prob_threshold,
) where {T}
    id = one(parent(first(polys)))
    return field_contains(field, [p // id for p in polys], prob_threshold)
end

# ------------------------------------------------------------------------------

function issubfield(
    F::RationalFunctionField{T},
    E::RationalFunctionField{T},
    prob_threshold,
) where {T}
    return all(field_contains(E, F.dennums, prob_threshold))
end

function fields_equal(
    F::RationalFunctionField{T},
    E::RationalFunctionField{T},
    prob_threshold,
) where {T}
    new_p = 1 - (1 - prob_threshold) / 2
    return issubfield(F, E, new_p) && issubfield(E, F, new_p)
end

# ------------------------------------------------------------------------------

"""
    beautiful_generators(rff::RationalFunctionField)

Given a field of rational functions `rff` returns a set of "simpler" and
standardized generators for `rff`.

Applies the following passes:
1. Filter constants,
2. Remove redundant generators.
"""
@timeit _to function beautiful_generators(
    rff::RationalFunctionField;
    discard_redundant = true,
    reversed_order = false,
    priority_variables = [],
)
    time_start = time_ns()
    fracs = dennums_to_fractions(rff.dennums)
    # Filter pass
    fracs = filter(!is_rational_func_const, fracs)
    fracs = unique(fracs)
    if isempty(fracs)
        @debug "The set of generators is empty"
        return fracs
    end
    # Remove redundant pass
    if discard_redundant
        fracs_priority = filter(f -> issubset(vars(f), priority_variables), fracs)
        fracs_rest = filter(f -> !(f in fracs_priority), fracs)
        sort!(fracs_priority, lt = rational_function_cmp)
        sort!(fracs_rest, lt = rational_function_cmp)
        fracs = vcat(fracs_priority, fracs_rest)
        @debug "The pool of fractions:\n$(join(map(repr, fracs), ",\n"))"
        if reversed_order
            non_redundant = collect(1:length(fracs))
            for i in length(fracs):-1:1
                func = fracs[i]
                if length(non_redundant) == 1
                    continue
                end
                result = field_contains_mod_p(
                    RationalFunctionField(fracs[setdiff(non_redundant, i)]),
                    [func],
                )
                @debug "Simplification: inclusion check $func $result"
                if result[1]
                    @debug "The function $func is discarded"
                    setdiff!(non_redundant, i)
                end
            end
        else
            non_redundant = Vector{Int}()
            push!(non_redundant, 1)
            for i in 2:length(fracs)
                func = fracs[i]
                result = field_contains_mod_p(
                    RationalFunctionField(fracs[non_redundant]),
                    [func],
                )
                @debug "Simplification: inclusion check $func $result"
                if !result[1]
                    @debug "The function $func is included in the set of generators"
                    push!(non_redundant, i)
                end
            end
        end
        @debug "Out of $(length(fracs)) simplified generators there are $(length(non_redundant)) non redundant"
        fracs = fracs[non_redundant]
    end
    sort!(fracs, lt = (f, g) -> rational_function_cmp(f, g))
    spring_cleaning_pass!(fracs)
    _runtime_logger[:id_beautifulization] += (time_ns() - time_start) / 1e9
    return fracs
end

function spring_cleaning_pass!(fracs)
    @assert all(is_rational_func_normalized, fracs)
    for i in 1:length(fracs)
        func = fracs[i]
        num, den = unpack_fraction(func)
        if is_constant(num)
            func = den // num
        end
        num, den = unpack_fraction(func)
        if leading_coefficient(num) < 0
            func = func * leading_coefficient(num)
        end
        num, den = unpack_fraction(func)
        if !isone(leading_coefficient(num))
            func = divexact(func, leading_coefficient(num))
        end
        num, den = unpack_fraction(func)
        if is_constant(den) && is_constant(Nemo.term(num, length(num)))
            func = (num - trailing_coefficient(num)) // one(num)
        end
        fracs[i] = func
    end
    fracs
end

# ------------------------------------------------------------------------------

"""
    groebner_basis_coeffs(rff; options...)

## Options

- `ordering`: GB ordering; must be one of the orderings exported by
  `ParamPunPam` or `Groebner`.
- `up_to_degree`: a tuple of integers, the degrees of numerator and denominator.
    The result is correct up to the requested degrees.
"""
@timeit _to function groebner_basis_coeffs(
    rff::RationalFunctionField;
    seed = 42,
    ordering = Groebner.InputOrdering(),
    up_to_degree = (typemax(Int), typemax(Int)),
    rational_interpolator = :VanDerHoevenLecerf,
)
    mqs = rff.mqs
    if are_generators_zero(mqs)
        return rff
    end
    gb, fracs, new_rff = nothing, nothing, nothing
    # Check if the basis is in cache
    if haskey(mqs.cached_groebner_bases, (ordering, up_to_degree))
        @debug "Cache hit with ($ordering, $up_to_degree)!"
        gb = mqs.cached_groebner_bases[ordering, up_to_degree]
        basis_coeffs = map(collect ∘ coefficients, gb)
        fracs = collect(mapreduce(Set, union!, basis_coeffs))
        return RationalFunctionField(fracs)
    end
    _runtime_logger[:id_calls_to_gb] += 1
    current_degrees = (2, 2)
    two_sided_inclusion = false
    while !two_sided_inclusion && all(current_degrees .<= up_to_degree)
        @debug "Computing GB with parameters up to degrees $(current_degrees)"
        runtime = @elapsed gb = ParamPunPam.paramgb(
            mqs,
            up_to_degree = current_degrees,
            ordering = ordering,
            rational_interpolator = rational_interpolator,
        )
        _runtime_logger[:id_npoints_degree] +=
            ParamPunPam._runtime_data[:npoints_degree_estimation]
        _runtime_logger[:id_npoints_interpolation] +=
            ParamPunPam._runtime_data[:npoints_interpolation]
        _runtime_logger[:id_groebner_time] += runtime
        @debug "Groebner basis computed in $runtime seconds"
        basis_coeffs = map(collect ∘ coefficients, gb)
        basis_coeffs_set = mapreduce(Set, union!, basis_coeffs)
        fracs = collect(basis_coeffs_set)
        @debug "Generators up to degrees $(current_degrees) are $fracs"
        @debug "Checking two-sided inclusion modulo a prime"
        time_start = time_ns()
        # Check inclusion: <simplified generators> in <original generators> 
        new_rff = RationalFunctionField(fracs)
        inclusion = issubfield_mod_p(new_rff, rff)
        two_sided_inclusion = two_sided_inclusion || inclusion
        # Check inclusion: <original generators> in <simplified generators>
        inclusion = issubfield_mod_p(rff, new_rff)
        runtime = (time_ns() - time_start) / 1e9
        _runtime_logger[:id_inclusion_check_mod_p] += runtime
        two_sided_inclusion = two_sided_inclusion && inclusion
        @debug "Inclusion checked in $(runtime) seconds. Result: $two_sided_inclusion"
        current_degrees = current_degrees .* 2
    end
    @debug "The coefficients of the Groebner basis are presented by $(length(fracs)) rational functions"
    new_rff.mqs.cached_groebner_bases[ordering, up_to_degree] = gb
    rff.mqs.cached_groebner_bases[ordering, up_to_degree] = gb
    return new_rff
end

"""
    generating_sets_fan(rff::RationalFunctionField, nbases)

Returns a set of Groebner bases for multiple different rankings of variables.

## Arguments

- `nbases`: How many bases to compute.
- Keyword `up_to_degree`: a tuple of integers, max. degrees of numerators and
  denominators. Result is correct up to the requested degrees.
"""
@timeit _to function generating_sets_fan(
    rff::RationalFunctionField{T},
    code::Integer;
    seed = 42,
    up_to_degree = (3, 3),
) where {T}
    time_start = time_ns()
    vars = gens(parent(rff.mqs))
    nbases = length(vars)
    @info "Computing $nbases Groebner bases for degrees $up_to_degree for block orderings"
    ordering_to_generators = Dict()
    if code == 0
        return ordering_to_generators
    end
    # The first basis in some ordering
    ord = InputOrdering()
    new_rff = groebner_basis_coeffs(rff, seed = seed, ordering = ord)
    cfs = beautiful_generators(new_rff)
    ordering_to_generators[ord] = cfs
    if isempty(cfs)
        return ordering_to_generators
    end
    if length(vars) == 1
        return ordering_to_generators
    end
    # NOTE: maybe hide the computation of multiple bases inside
    # RationalFunctionField
    gb_rff = RationalFunctionField(cfs)
    if code >= 1
        for i in 1:nbases
            vars_shuffled = circshift(vars, i)
            n = length(vars_shuffled)
            # n1, n2 = div(n, 2), n - div(n, 2)
            n1, n2 = n - 1, 1
            ord = DegRevLex(vars_shuffled[1:n1]) * DegRevLex(vars_shuffled[(n1 + 1):end])
            @debug "Computing GB for ordering $ord"
            new_rff = groebner_basis_coeffs(
                gb_rff,
                seed = seed,
                ordering = ord,
                up_to_degree = up_to_degree,
            )
            cfs = beautiful_generators(new_rff, discard_redundant = false)
            ordering_to_generators[ord] = cfs
        end
    end
    if code >= 2
        for _ in 1:nbases
            vars_shuffled = shuffle(vars)
            n = length(vars_shuffled)
            n1, n2 = max(n - 2, 1), min(2, length(vars) - 1)
            ord = DegRevLex(vars_shuffled[1:n1]) * DegRevLex(vars_shuffled[(n1 + 1):end])
            @debug "Computing GB for ordering $ord"
            new_rff = groebner_basis_coeffs(
                gb_rff,
                seed = seed,
                ordering = ord,
                up_to_degree = up_to_degree,
            )
            cfs = beautiful_generators(new_rff, discard_redundant = false)
            ordering_to_generators[ord] = cfs
        end
    end
    if code >= 3
        for _ in 1:nbases
            vars_shuffled = shuffle(vars)
            n = length(vars_shuffled)
            n1 = div(n, 2)
            n2 = n - n1
            ord = DegRevLex(vars_shuffled[1:n1]) * DegRevLex(vars_shuffled[(n1 + 1):end])
            @debug "Computing GB for ordering $ord"
            new_rff = groebner_basis_coeffs(
                gb_rff,
                seed = seed,
                ordering = ord,
                up_to_degree = up_to_degree,
            )
            cfs = beautiful_generators(new_rff, discard_redundant = false)
            ordering_to_generators[ord] = cfs
        end
    end
    _runtime_logger[:id_gbfan_time] = (time_ns() - time_start) / 1e9
    @info "Computed Groebner bases in $((time_ns() - time_start) / 1e9) seconds"
    return ordering_to_generators
end

function monomial_generators_up_to_degree(
    rff::RationalFunctionField{T},
    up_to_degree;
    seed = 42,
    strategy = :monte_carlo,
) where {T}
    @assert strategy in (:monte_carlo,)
    relations = linear_relations_between_normal_forms(
        beautiful_generators(rff),
        up_to_degree,
        seed = seed,
    )
    return relations
end

"""
    simplified_generating_set(rff; prob_threshold = 0.99, seed = 42)

Returns a simplified set of generators for `rff`. 
Result is correct (in the Monte-Carlo sense) with probability at least `prob_threshold`.
"""
@timeit _to function simplified_generating_set(
    rff::RationalFunctionField;
    prob_threshold = 0.99,
    seed = 42,
    simplify = :standard,
    check_variables = false, # almost always slows down and thus turned off
    rational_interpolator = :VanDerHoevenLecerf,
    priority_variables = [],
)
    @info "Simplifying generating set. Simplification level: $simplify"
    _runtime_logger[:id_groebner_time] = 0.0
    _runtime_logger[:id_calls_to_gb] = 0
    _runtime_logger[:id_inclusion_check_mod_p] = 0.0
    _runtime_logger[:id_inclusion_check] = 0.0
    _runtime_logger[:id_gbfan_time] = 0.0
    _runtime_logger[:id_normalforms_time] = 0.0
    _runtime_logger[:id_ranking] = 0

    # Checking membership of particular variables and adding them to the field
    if check_variables
        vars = gens(poly_ring(rff))
        containment = field_contains(rff, vars, (1.0 + prob_threshold) / 2)
        prob_threshold = (1.0 + prob_threshold) / 2
        if all(containment)
            return [v // one(poly_ring(rff)) for v in vars]
        end
        field_gens = rff.dennums
        for (v, is_contained) in zip(vars, containment)
            if is_contained
                push!(field_gens, [one(poly_ring(rff)), v])
            end
        end
        rff = RationalFunctionField(field_gens)
    end

    normalforms_degree = 2
    gbfan_simplification_code = 1
    if simplify === :standard
        # pass
    elseif simplify === :weak
        normalforms_degree = 2
        gbfan_simplification_code = 0
    elseif simplify === :strong
        normalforms_degree = 3
        gbfan_simplification_code = 3
    end

    # Compute the first GB in some ordering
    new_rff = groebner_basis_coeffs(
        rff,
        seed = seed,
        rational_interpolator = rational_interpolator,
    )
    new_fracs = beautiful_generators(new_rff)
    if isempty(new_fracs)
        return new_fracs
    end

    # Compute some normal forms
    start_time = time_ns()
    rff_generators = monomial_generators_up_to_degree(
        new_rff,
        normalforms_degree;
        seed = seed,
        strategy = :monte_carlo,
    )
    append!(new_fracs, rff_generators)
    _runtime_logger[:id_normalforms_time] = (time_ns() - start_time) / 1e9

    # Compute some GBs
    fan = generating_sets_fan(new_rff, gbfan_simplification_code; seed = seed)
    for (ord, rff_gens) in fan
        append!(new_fracs, rff_gens)
    end
    new_fracs_unique = unique(new_fracs)
    @debug """
Final cleaning and simplification of generators. 
Out of $(length(new_fracs)) fractions $(length(new_fracs_unique)) are syntactically unique."""
    runtime = @elapsed new_fracs = beautiful_generators(
        RationalFunctionField(new_fracs_unique),
        priority_variables = priority_variables,
    )
    @debug "Checking inclusion with probability $prob_threshold"
    runtime =
        @elapsed result = issubfield(rff, RationalFunctionField(new_fracs), prob_threshold)
    _runtime_logger[:id_inclusion_check] = runtime
    if !result
        @warn "Field membership check failed. Error will follow."
        throw("The new subfield generators are not correct.")
    end
    @info "Inclusion checked with probability $prob_threshold in $(_runtime_logger[:id_inclusion_check]) seconds"
    @debug "Out of $(length(rff.mqs.nums_qq)) initial generators there are $(length(new_fracs)) independent"
    ranking = generating_set_rank(new_fracs)
    _runtime_logger[:id_ranking] = ranking
    @debug "The ranking of the new set of generators is $ranking"
    return new_fracs
end

# ------------------------------------------------------------------------------
