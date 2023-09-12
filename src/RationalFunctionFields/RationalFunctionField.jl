
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

    function RationalFunctionField(polys::Vector{T}) where {T}
        RationalFunctionField(polys .// one(parent(first(polys))))
    end
    function RationalFunctionField(fractions::Vector{Generic.Frac{T}}) where {T}
        RationalFunctionField(fractions_to_dennums(fractions))
    end
    function RationalFunctionField(dennums::Vector{Vector{T}}) where {T}
        @assert !isempty(dennums)
        new{T}(dennums, IdealMQS(dennums))
    end
end

# ------------------------------------------------------------------------------

function poly_ring(F::RationalFunctionField)
    return parent(first(first(F.dennums)))
end

# ------------------------------------------------------------------------------

"""
    check_field_membership_mod_p(generators, rat_funcs)

Checks whether given rational functions belong to a given field of rational
functions over integers modulo a prime

Inputs:
- `generators` - a list of lists of polynomials. Each of the lists, say, `[f1, ..., fn]`,
  defines generators `f2/f1, ..., fn/f1`. Let ``F`` be the field generated by all of them.
- `rat_funcs` - list of rational functions

Output:
- a list `L[i]` of bools of length `length(rat_funcs)` such that `L[i]` is true iff
   the i-th function belongs to ``F``
"""
function check_field_membership_mod_p(generators, rat_funcs)
    if isempty(generators)
        # TODO: RationalFunctionField of empty set is an error currently
        return fill(false, length(rat_funcs))
    end
    if isempty(rat_funcs)
        return Bool[]
    end

    check_field_membership_mod_p!(
        RationalFunctionField(generators),
        RationalFunctionField(rat_funcs),
    )
end

# ------------------------------------------------------------------------------

"""
    field_contains(field, ratfuncs, p)

Checks whether given rational function field `field` contains given rational
functions `ratfuncs` (represented as a list of lists). The result is correct with
probability at least `p`

Inputs:
- `field` - a rational function field
- `ratfuncs` - a list of lists of polynomials. Each of the lists, say, `[f1, ..., fn]`,
  defines generators `f2/f1, ..., fn/f1`.
- `p` real number from (0, 1)

Output:
- a list `L[i]` of bools of length `length(rat_funcs)` such that `L[i]` is true iff
   the i-th function belongs to `field`
"""
function field_contains(
    field::RationalFunctionField{T},
    ratfuncs::Vector{Vector{T}},
    p,
) where {T}
    if isempty(ratfuncs)
        return Bool[]
    end
    @debug "Finding pivot polynomials"
    pivots = map(plist -> plist[findmin(map(total_degree, plist))[2]], ratfuncs)
    @debug "\tDegrees are $(map(total_degree, pivots))"

    @debug "Estimating the sampling bound"
    # uses Theorem 3.3 from https://arxiv.org/pdf/2111.00991.pdf
    # the comments below use the notation from the theorem
    ring = parent(first(first(ratfuncs)))
    den_lcm = lcm(field.mqs.den_lcm_orig, foldl(lcm, pivots))
    # this is deg(g) + 1
    degree = total_degree(den_lcm) + 1
    # computing maximum of deg(f) for different f's to be tested
    for (i, plist) in enumerate(ratfuncs)
        extra_degree = total_degree(den_lcm) - total_degree(pivots[i])
        degree = max(degree, extra_degree + maximum(total_degree, plist))
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
        (length(ratfuncs)) *
        ceil(1 / (1 - p)),
    )
    @debug "\tSampling from $(-sampling_bound) to $(sampling_bound)"

    mqs = field.mqs
    param_ring = ParamPunPam.parent_params(mqs)
    point = map(v -> Nemo.QQ(rand((-sampling_bound):sampling_bound)), gens(param_ring))
    mqs_specialized = specialize(mqs, point)
    @debug "Computing Groebner basis ($(length(mqs_specialized)) equations)"
    mqs_ratfuncs = specialize(IdealMQS(ratfuncs), point; saturated = false)
    @assert parent(first(mqs_specialized)) == parent(first(mqs_ratfuncs))
    gb = groebner(mqs_specialized)
    result = map(iszero, normalform(gb, mqs_ratfuncs))
    return result
end

function field_contains(
    field::RationalFunctionField{T},
    ratfuncs::Vector{Generic.Frac{T}},
    p,
) where {T}
    return field_contains(field, fractions_to_dennums(ratfuncs), p)
end

function field_contains(field::RationalFunctionField{T}, polys::Vector{T}, p) where {T}
    id = one(parent(first(polys)))
    return field_contains(field, [[id, p] for p in polys], p)
end

# ------------------------------------------------------------------------------

function issubfield(F::RationalFunctionField{T}, E::RationalFunctionField{T}, p) where {T}
    return all(field_contains(E, F.dennums, p))
end

function fields_equal(F::RationalFunctionField{T}, E::RationalFunctionField{T}, p) where {T}
    new_p = 1 - (1 - p) / 2
    return issubfield(F, E, new_p) && issubfield(E, F, new_p)
end

# ------------------------------------------------------------------------------

function check_field_membership_mod_p!(
    generators::RationalFunctionField{T},
    tobereduced::RationalFunctionField{T},
) where {T}
    mqs_generators = generators.mqs
    mqs_tobereduced = tobereduced.mqs
    ff = Nemo.GF(2^31 - 1)
    reduce_mod_p!(mqs_generators, ff)
    reduce_mod_p!(mqs_tobereduced, ff)
    param_ring = ParamPunPam.parent_params(mqs_generators)
    point = ParamPunPam.distinct_nonzero_points(ff, nvars(param_ring))
    gens_specialized = ParamPunPam.specialize_mod_p(mqs_generators, point)
    polys_specialized =
        ParamPunPam.specialize_mod_p(mqs_tobereduced, point, saturated = false)
    @assert parent(first(gens_specialized)) == parent(first(polys_specialized))
    gb = groebner(gens_specialized)
    nf = normalform(gb, polys_specialized)
    result = map(iszero, nf)
    return result
end

# ------------------------------------------------------------------------------

"""
    beautifuly_generators(rff::RationalFunctionField)

Given a field of rational functions `rff` returns a set of "simpler" and
standardized generators for `rff`.

Applies the following passes:
1. Filter constants,
2. Remove redundant generators.
"""
function beautifuly_generators(
    rff::RationalFunctionField;
    discard_redundant = true,
    reversed_order = false,
)
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
        sort!(fracs, lt = rational_function_cmp)
        @debug "The pool of fractions:\n$(join(map(repr, fracs), ",\n"))"
        if reversed_order
            non_redundant = collect(1:length(fracs))
            for i in length(fracs):-1:1
                func = fracs[i]
                if length(non_redundant) == 1
                    continue
                end
                result =
                    check_field_membership_mod_p(fracs[setdiff(non_redundant, i)], [func])
                @debug "Simplification: inclusion check" func result
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
                result = check_field_membership_mod_p(fracs[non_redundant], [func])
                @debug "Simplification: inclusion check" func result
                if !result[1]
                    @debug "The function $func is included in the set of generators"
                    push!(non_redundant, i)
                end
            end
        end
        @debug "Out of $(length(fracs)) simplified generators there are $(length(non_redundant)) non redundant"
        fracs = fracs[non_redundant]
    end
    sort!(fracs, lt = rational_function_cmp)
    spring_cleaning_pass!(fracs)
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
function groebner_basis_coeffs(
    rff::RationalFunctionField;
    seed = 42,
    ordering = Groebner.InputOrdering(),
    up_to_degree = (typemax(Int), typemax(Int)),
    use_homogenization = false,
    rational_interpolator = :VanDerHoevenLecerf,
)
    mqs = rff.mqs
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
        @info "Groebner basis computed in $runtime seconds"
        basis_coeffs = map(collect ∘ coefficients, gb)
        basis_coeffs_set = mapreduce(Set, union!, basis_coeffs)
        fracs = collect(basis_coeffs_set)
        @debug "Generators up to degrees $(current_degrees) are" fracs
        @info "Checking two-sided inclusion modulo a prime"
        time_start = time_ns()
        # Check inclusion: <simplified generators> in <original generators> 
        new_rff = RationalFunctionField(fracs)
        inclusion = check_field_membership_mod_p!(rff, new_rff)
        two_sided_inclusion = two_sided_inclusion || all(inclusion)
        # Check inclusion: <original generators> in <simplified generators>
        inclusion = check_field_membership_mod_p!(new_rff, rff)
        runtime = (time_ns() - time_start) / 1e9
        _runtime_logger[:id_inclusion_check_mod_p] += runtime
        @info "Inclusion checked in $(runtime) seconds. Result: $two_sided_inclusion"
        two_sided_inclusion = two_sided_inclusion && all(inclusion)
        current_degrees = current_degrees .* 2
    end
    @info "The coefficients of the Groebner basis are presented by $(length(fracs)) rational functions"
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
function generating_sets_fan(
    rff::RationalFunctionField{T},
    code::Integer;
    seed = 42,
    up_to_degree = (3, 3),
) where {T}
    time_start = time_ns()
    vars = gens(parent(rff.mqs))
    nbases = length(vars)
    @info "Computing $nbases Groebner bases for block orderings. Simplification code is $code"
    ordering_to_generators = Dict()
    # The first basis in some ordering
    ord = InputOrdering()
    new_rff = groebner_basis_coeffs(rff, seed = seed, ordering = ord)
    cfs = beautifuly_generators(new_rff)
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
            @info "Computing GB for ordering" ord
            new_rff = groebner_basis_coeffs(
                gb_rff,
                seed = seed,
                ordering = ord,
                up_to_degree = up_to_degree,
                use_homogenization = true,
            )
            cfs = beautifuly_generators(new_rff, discard_redundant = false)
            ordering_to_generators[ord] = cfs
        end
    end
    if code >= 2
        for _ in 1:nbases
            vars_shuffled = shuffle(vars)
            n = length(vars_shuffled)
            n1, n2 = max(n - 2, 1), min(2, length(vars) - 1)
            ord = DegRevLex(vars_shuffled[1:n1]) * DegRevLex(vars_shuffled[(n1 + 1):end])
            @info "Computing GB for ordering" ord
            new_rff = groebner_basis_coeffs(
                gb_rff,
                seed = seed,
                ordering = ord,
                up_to_degree = up_to_degree,
                use_homogenization = true,
            )
            cfs = beautifuly_generators(new_rff, discard_redundant = false)
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
            @info "Computing GB for ordering" ord
            new_rff = groebner_basis_coeffs(
                gb_rff,
                seed = seed,
                ordering = ord,
                up_to_degree = up_to_degree,
                use_homogenization = true,
            )
            cfs = beautifuly_generators(new_rff, discard_redundant = false)
            ordering_to_generators[ord] = cfs
        end
    end
    _runtime_logger[:id_gbfan_time] = (time_ns() - time_start) / 1e9
    return ordering_to_generators
end

function monomial_generators_up_to_degree(
    rff::RationalFunctionField{T},
    up_to_degree;
    seed = 42,
    strategy = :monte_carlo,
) where {T}
    @assert strategy in (:monte_carlo,)
    relations = linear_relations_between_normal_forms_mod_p(
        beautifuly_generators(rff),
        up_to_degree,
        seed = seed,
    )
    return relations
end

"""
    simplified_generating_set(rff; p = 0.99, seed = 42)

Returns a simplified set of generators for `rff`. 
Result is correct (in Monte-Carlo sense) with probability at least `p`.
"""
function simplified_generating_set(
    rff::RationalFunctionField;
    p = 0.99,
    seed = 42,
    strategy = (:gb,),
    check_variables = false, # almost always slows down and thus turned off
    rational_interpolator = :VanDerHoevenLecerf,
)
    @info "Simplifying identifiable functions"
    _runtime_logger[:id_groebner_time] = 0.0
    _runtime_logger[:id_calls_to_gb] = 0
    _runtime_logger[:id_inclusion_check_mod_p] = 0.0
    _runtime_logger[:id_inclusion_check] = 0.0
    _runtime_logger[:id_beautifulization] = 0.0
    _runtime_logger[:id_gbfan_time] = 0.0
    _runtime_logger[:id_normalforms_time] = 0.0
    _runtime_logger[:id_ranking] = 0

    # Checking identifiability of particular variables and adding them to the field
    if check_variables
        vars = gens(poly_ring(rff))
        containment = field_contains(rff, vars, (1.0 + p) / 2)
        p = (1.0 + p) / 2
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

    # Compute the first GB in some ordering
    new_rff = groebner_basis_coeffs(
        rff,
        seed = seed,
        rational_interpolator = rational_interpolator,
    )
    new_fracs = beautifuly_generators(new_rff)
    if isempty(new_fracs)
        return new_fracs
    end
    # If a set of GBs is needed
    if first(strategy) === :gbfan
        @assert length(strategy) == 2
        _, nbases = strategy
        fan = generating_sets_fan(new_rff, nbases; seed = seed)
        for (ord, generators) in fan
            append!(new_fracs, generators)
        end
    end
    # If normal forms are needed
    if first(strategy) === :normalforms
        @assert length(strategy) == 2
        _, up_to_degree = strategy
        generators = monomial_generators_up_to_degree(
            new_rff,
            up_to_degree;
            seed = seed,
            strategy = :monte_carlo,
        )
        append!(new_fracs, generators)
    end
    # Something in the middle
    if first(strategy) === :hybrid
        @assert 1 <= length(strategy) <= 2
        if length(strategy) == 1
            nbases = 10
        else
            nbases = strategy[2]
        end

        # Compute some normal forms
        up_to_degree = 3
        generators = monomial_generators_up_to_degree(
            new_rff,
            up_to_degree;
            seed = seed,
            strategy = :monte_carlo,
        )
        append!(new_fracs, generators)

        # Compute some GBs
        fan = generating_sets_fan(new_rff, nbases; seed = seed)
        for (ord, generators) in fan
            append!(new_fracs, generators)
        end
    end
    new_fracs_unique = unique(new_fracs)
    @info """
    Final cleaning and simplification of generators. 
    Out of $(length(new_fracs)) fractions $(length(new_fracs_unique)) are syntactically unique."""
    runtime =
        @elapsed new_fracs = beautifuly_generators(RationalFunctionField(new_fracs_unique))
    _runtime_logger[:id_beautifulization] += runtime
    @info "Checking inclusion with probability $p"
    runtime = @elapsed result = issubfield(rff, RationalFunctionField(new_fracs), p)
    _runtime_logger[:id_inclusion_check] = runtime
    @info "Inclusion checked in $(_runtime_logger[:id_inclusion_check]) seconds. Result: $result"
    if !result
        @warn "Field membership check failed. Error will follow."
        throw("The new subfield generators are not correct.")
    end
    @info "Out of $(length(rff.mqs.nums_qq)) initial generators there are $(length(new_fracs)) indepdendent"
    ranking = generating_set_rank(new_fracs)
    _runtime_logger[:id_ranking] = ranking
    @info "The ranking of the new set of generators is $ranking"
    return new_fracs
end

# ------------------------------------------------------------------------------
