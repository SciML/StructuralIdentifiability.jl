"""
    check_constructive_field_membership(tagged_mqs, tag_relations, tagged_num, tagged_den)

Reduces `tagged_num // tagged_den` modulo the given `tagged_mqs` ideal.
Cancels out `tag_relations`.

## Input

- `tagged_mqs`: vector of generators in K(T)[x].
- `tag_relations`: vector of relations in K[T].
- `tagged_num`, `tagged_den`: numerator and denominator to be reduced, both
  elements of `K(T)[x]`.

## Output

Return a tuple (`membership`, `remainder`).

- `membership`: `true` if `tagged_num // tagged_den` belongs to the field
    associated with `tagged_mqs`, `false` otherwise.
- `remainder`: expression of `tagged_num // tagged_den` in terms of `T`.
    NOTE: expression is canonical provided that `tag_relations` are specified.
"""
function check_constructive_field_membership(
    tagged_mqs,
    tag_relations,
    tagged_num,
    tagged_den,
)
    @assert !isempty(tagged_mqs)
    @assert !iszero(tagged_den)
    ring_of_tags = base_ring(base_ring(parent(first(tagged_mqs))))
    if !isempty(tag_relations)
        @assert ring_of_tags == parent(first(tag_relations))
    end
    # Compute the remainders modulo the MQS.
    # These belong to K(T)[x].
    @debug """
    Reducing $tagged_num, $tagged_den"""
    _, num_rem = divrem(tagged_num, tagged_mqs)
    _, den_rem = divrem(tagged_den, tagged_mqs)
    @debug """
    Normal forms modulo MQS: 
    Num: $(num_rem)
    Den: $(den_rem)"""

    # reducing the coefficients modulo GB
    num_rem = normalize_coefficients(num_rem, tag_relations)
    den_rem = normalize_coefficients(den_rem, tag_relations)
    if iszero(den_rem)
        @debug """
        The element $tagged_num // $tagged_den is not in the sub-field
        Normal form, numerator: $num_rem
        Normal form, denominator: $den_rem
        """
        return false, zero(ring_of_tags) // one(ring_of_tags)
    end

    m = first(monomials(den_rem))
    c_den = first(coefficients(den_rem))
    c_num = coeff(num_rem, m)
    fraction_candidate = c_num // c_den
    rem_rem = num_rem - fraction_candidate * den_rem
    rem_rem = normalize_coefficients(rem_rem, tag_relations)

    # If norm_form(Num) // norm_form(Den) does not belongs to K(T), then
    # the fraction does not belong to the field
    if !iszero(rem_rem)
        @debug """
        The element $tagged_num // $tagged_den is not in the sub-field
        Normal form, numerator: $num_rem
        Normal form, denominator: $den_rem
        Remainder: $rem_rem
        """
        return false, zero(ring_of_tags) // one(ring_of_tags)
    end
    # Now we know that the remainders belong to K(T). 
    # To obtain a canonical representation We need to cancel out the algebraic
    # relations in K[T] between the tags
    num_tags = numerator(fraction_candidate)
    den_tags = denominator(fraction_candidate)
    _, num_tags_factored = divrem(num_tags, tag_relations)
    _, den_tags_factored = divrem(den_tags, tag_relations)
    @debug """
    Before factoring out relations:
    Num: $(num_tags)
    Den: $(den_tags)
    After factoring out relations:
    Num: $(num_tags_factored)
    Den: $(den_tags_factored)
    """
    rem_canon = num_tags_factored // den_tags_factored
    return true, rem_canon
end

"""
    check_constructive_field_membership(rff, to_be_reduced)

Returns the unique expression of each fraction in `to_be_reduced` in terms of
the elements of the given rational function field `rff`.

Follows the vein of Algorithm 1.17 from https://doi.org/10.1006/jsco.1998.0246

## Input

- `rff`: a subfield of the field of rational functions.
- `to_be_reduced`: an array of fractions to be reduced.
- `tag_names` (optional): a vector of strings to be used as tags.

## Output

Returns a 4-element tuple 
(`memberships`, `remainders`, `relations_between_tags`, `tag_to_gen`).

- `memberships`: is `true` whenever `to_be_reduced[i]` belongs to the field.
- `remainders`: is the unique expression of `to_be_reduced[i]` in terms of the
    generators of `rff`. If the membership was false, nothing undefined behaviour
- `relations_between_tags`: is the list of algebraic relations between the
    generators of `rff`.
- `tag_to_gen`: is a dictionary mapping each new tag to the corresponding
    generator of `rff`.

"""
function check_constructive_field_membership(
    rff::RationalFunctionField{T},
    to_be_reduced::Vector{Generic.FracFieldElem{T}};
    tag_names = Vector{String}(),
) where {T}
    @assert !isempty(to_be_reduced)
    fracs_gen = generators(rff)
    @assert parent(first(fracs_gen)) == parent(first(to_be_reduced))
 
    algebraicity = check_algebraicity(rff, to_be_reduced, 0.99)
    to_be_reduced = to_be_reduced[algebraicity]
    # A tag is assigned for each of the the generators of the given rational
    # function field. Then, the MQS ideal is constructed in the following way:
    #   < N_i(x) - T_i * D_i(x), Q * t - 1 >  in  K(T)[x][t]
    # for each N_i / D_i in the generators.
    #
    # Let the fraction to be reduced be Num // Den, A is a formal parameter.
    # We construct Num - A * Den in K[x][A] to later compute the normal form
    # of it with respect to the MQS ideal of generators in K(T)[vars][A]. 
    #
    # Then, let A = norm_form(Num) // norm_form(Den). Note that A lives in
    # K(T)[x].
    #
    # If A belongs to K(T), then Num // Den belongs to the the given field.
    # Otherwise, A belongs to K(T)[x], and Num // Den does not belong to the
    # given field.
    x_ring = poly_ring(rff)
    K = base_ring(x_ring)
    orig_strings = map(string, gens(x_ring))
    # add extra dummy tags for the elements of the transcendence basis
    # the idea is that they will not appear in the final expressions anyway
    tag_strings = if !isempty(tag_names)
        @assert length(fracs_gen) == length(tag_names)
        vcat(tag_names, gen_tag_names(length(rff.trbasis_over), "Tag"))
    else
        gen_tag_names(length(fracs_gen) + length(rff.trbasis_over), "Tag")
    end
    sat_string = gen_tag_name("Sat")
    @debug """
Tags:
$(join(map(x -> string(x[1]) * " -> " * string(x[2]),  zip(fracs_gen, tag_strings)), "\t\n"))
Saturation tag:
$sat_string
"""
    extended_gens = vcat(fracs_gen, rff.trbasis_over)
    poly_ring_tag, vars_tag =
        polynomial_ring(K, vcat(sat_string, orig_strings, tag_strings))
    sat_var = vars_tag[1]
    orig_vars = vars_tag[2:(nvars(x_ring) + 1)]
    tag_vars = vars_tag[(nvars(x_ring) + 2):end]
    # Construct generators of the tagged MQS ideal.
    tagged_mqs = Vector{elem_type(poly_ring_tag)}(undef, length(extended_gens) + 1)
    Q = one(poly_ring_tag)
    for i in 1:length(extended_gens)
        num, den = unpack_fraction(extended_gens[i])
        num_tag = parent_ring_change(num, poly_ring_tag)
        den_tag = parent_ring_change(den, poly_ring_tag)
        Q = lcm(Q, squarefree_part(den_tag))
        tagged_poly_mqs = num_tag - tag_vars[i] * den_tag
        tagged_mqs[i] = tagged_poly_mqs
    end
    tagged_mqs[end] = Q * sat_var - 1
    # Compute the basis of the MQS in K[T][x][t] such that T < x < t.
    #
    # NOTE: we compute the basis in K[T][x][t], not in K(T)[x][t].
    # This way, we obtain two pieces of information at once:
    # - the kernel of the map from T to x (the algebraic relations of the tags).
    # - the GB of the MQS in K(T)[x][t].
    # ord = Lex()
    ord = DegRevLex([sat_var]) * DegRevLex(orig_vars) * DegRevLex(tag_vars)
    @debug """
    Tagged MQS ideal:
    $tagged_mqs
    Monom ordering:
    $(ord)"""
    tagged_mqs_gb = groebner(tagged_mqs, ordering = ord, homogenize = :no)
    # Relations between tags in K[T]
    relations_between_tags = filter(
        poly -> isempty(intersect(vars(poly), vcat(sat_var, orig_vars))),
        tagged_mqs_gb,
    )
    # The basis in K[T][x]
    tagged_mqs_gb = setdiff(tagged_mqs_gb, relations_between_tags)
    tagged_mqs_gb = filter(poly -> isempty(intersect(vars(poly), [sat_var])), tagged_mqs_gb)
    @debug """
    Tagged MQS GB:
    $tagged_mqs_gb
    Relations between tags:
    $relations_between_tags
    """
    # Reduce the fractions with respect to the MQS ideal.
    #
    # NOTE: reduction actually happens in K(T)[x]. So we map polynomials to the
    # parametric ring K(T)[x].
    ring_of_tags, tags = polynomial_ring(K, tag_strings)
    if !isempty(intersect(tag_strings, orig_strings))
        @warn """
    There is an intersection between the names of the tag variables and the original variables.
    Tags: $tag_strings
    Original vars: $orig_strings"""
    end
    parametric_ring, _ = polynomial_ring(
        fraction_field(ring_of_tags),
        orig_strings,
        internal_ordering = :degrevlex,
    )
    relations_between_tags =
        map(poly -> parent_ring_change(poly, ring_of_tags), relations_between_tags)
    param_var_mapping = merge(
        Dict(gens(poly_ring_tag)[2:(nvars(x_ring) + 1)] .=> gens(parametric_ring)),
        Dict(gens(poly_ring_tag)[(nvars(x_ring) + 2):end] .=> (gens(ring_of_tags)) .* one(parametric_ring)),
    )
    @debug """
    Variable mapping:
    $param_var_mapping
    Parametric ring:
    $parametric_ring
    """
    tagged_mqs_gb_param = map(
        poly -> eval_at_dict(poly, param_var_mapping),
        tagged_mqs_gb,
    )
    tagged_mqs_gb_param = map(f -> divexact(f, leading_coefficient(f)), tagged_mqs_gb_param)
    @debug "Tagged parametric mqs: $tagged_mqs_gb_param"
    # Reduce each fraction
    var_mapping = Dict(gens(x_ring) .=> gens(parametric_ring))
    memberships = Vector{Bool}(undef, length(to_be_reduced))
    remainders = Vector{Generic.FracFieldElem{T}}(undef, length(to_be_reduced))
    for i in 1:length(to_be_reduced)
        frac = to_be_reduced[i]
        num = eval_at_dict(numerator(frac), var_mapping)
        den = eval_at_dict(denominator(frac), var_mapping)
        membership, remainder = check_constructive_field_membership(
            tagged_mqs_gb_param,
            relations_between_tags,
            num,
            den,
        )
        memberships[i] = membership
        remainders[i] = remainder
    end

    # Cleaning up the "imaginary" tags corresponding to the transcendence basis
    short_ring_of_tags, tags = polynomial_ring(K, tag_strings[1:length(fracs_gen)])
    remainders = map(
        f ->
            parent_ring_change(numerator(f), short_ring_of_tags) //
            parent_ring_change(denominator(f), short_ring_of_tags),
        remainders,
    )
    new_remainders = Vector{Generic.FracFieldElem{T}}(undef, length(algebraicity))
    i = 1
    for j in 1:length(algebraicity)
        if !algebraicity[j]
            new_remainders[j] = one(short_ring_of_tags) // one(short_ring_of_tags)
        else
            new_remainders[j] = remainders[i]
            i += 1
        end
    end
    relations_between_tags =
        map(p -> parent_ring_change(p, short_ring_of_tags), relations_between_tags)
    tag_to_gen = Dict(tags[i] => fracs_gen[i] for i in 1:length(fracs_gen))

    return merge_results(algebraicity, memberships),
    new_remainders,
    relations_between_tags,
    tag_to_gen
end
