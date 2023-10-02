# Reparametrize!

# Maps the variables of the given polynomial according to `var_mapping`.
# The output polynomial lives in the given `new_ring`.
function crude_parent_ring_change(poly, new_ring, var_mapping)
    new_poly = zero(new_ring)
    for (i, term) in enumerate(terms(poly))
        new_term = one(new_ring) * coeff(poly, i)
        for var in vars(term)
            exp = degree(term, var)
            exp == 0 && continue
            new_var = var_mapping[var]
            new_term *= new_var^exp
        end
        new_poly += new_term
    end
    return new_poly
end

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
    # Compute the remainders module the MQS.
    # These belong to K(T)[x].
    @info """
    Reducing $tagged_num, $tagged_den"""
    _, num_rem = divrem(tagged_num, tagged_mqs)
    _, den_rem = divrem(tagged_den, tagged_mqs)
    @info """
    Normal forms modulo MQS: 
    Num: $(num_rem)
    Den: $(den_rem)"""
    common_factor = gcd(num_rem, den_rem)
    num_rem = divexact(num_rem, common_factor)
    den_rem = divexact(den_rem, common_factor)
    # If norm_form(Num) // norm_form(Den) does not belongs to K(T), then
    # the fraction does not belong to the field
    if iszero(den_rem)
        @warn """
        The element $tagged_num // $tagged_den is not in the sub-field
        Normal form, numerator: $num_rem
        Normal form, denominator: $den_rem
        Common factor: $(common_factor)
        """
        return false, zero(ring_of_tags) // one(ring_of_tags)
    end
    if total_degree(num_rem) > 0 || total_degree(den_rem) > 0
        @warn """
        The element $tagged_num // $tagged_den is not in the sub-field
        Normal form, numerator: $num_rem
        Normal form, denominator: $den_rem
        Common factor: $(common_factor)
        """
        return false, zero(ring_of_tags) // one(ring_of_tags)
    end
    # Now we know that the remainders belong to K(T). 
    # To obtain a canonical representation We need to cancel out the algebraic
    # relations in K[T] between the tags
    num = !iszero(num_rem) ? coeff(num_rem, 1) : zero(ring_of_tags) // one(ring_of_tags)
    den = coeff(den_rem, 1)
    num_num = numerator(num)
    num_den = numerator(den)
    _, num_num_factored = divrem(num_num, tag_relations)
    _, num_den_factored = divrem(num_den, tag_relations)
    num_factored = num_num_factored // denominator(num)
    den_factored = num_den_factored // denominator(den)
    @info """
    After factoring out relations:
    Num: $(num_factored)
    Den: $(den_factored)
    """
    if iszero(den_factored)
        return false, zero(ring_of_tags) // one(ring_of_tags)
    end
    rem_canon = num_factored // den_factored
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
    generators of `rff`.
- `relations_between_tags`: is the list of algebraic relations between the
    generators of `rff`.
- `tag_to_gen`: is a dictionary mapping each new tag to the corresponding
    generator of `rff`.

"""
function check_constructive_field_membership(
    rff::RationalFunctionField{T},
    to_be_reduced::Vector{Generic.Frac{T}};
    tag_names = Vector{String}(),
) where {T}
    @assert !isempty(to_be_reduced)
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
    fracs_gen = dennums_to_fractions(rff.dennums)
    @assert parent(first(fracs_gen)) == parent(first(to_be_reduced))
    poly_ring = base_ring(parent(first(to_be_reduced)))
    K = base_ring(poly_ring)
    orig_strings = map(string, gens(poly_ring))
    tag_strings = if !isempty(tag_names)
        @assert length(fracs_gen) == length(tag_names)
        tag_names
    else
        gen_tag_names(length(fracs_gen), "Tag")
    end
    sat_string = gen_tag_name("Sat")
    @info """
    Tags:
    $(join(map(x -> string(x[1]) * " -> " * string(x[2]),  zip(fracs_gen, tag_strings)), "\t\n"))
    Saturation tag:
    $sat_string
    """
    poly_ring_tag, vars_tag = PolynomialRing(K, vcat(sat_string, orig_strings, tag_strings))
    sat_var = vars_tag[1]
    orig_vars = vars_tag[2:(nvars(poly_ring) + 1)]
    tag_vars = vars_tag[(nvars(poly_ring) + 2):end]
    # Construct generators of the tagged MQS ideal.
    tagged_mqs = Vector{elem_type(poly_ring_tag)}(undef, length(fracs_gen) + 1)
    Q = one(poly_ring_tag)
    for i in 1:length(fracs_gen)
        num, den = unpack_fraction(fracs_gen[i])
        num_tag = parent_ring_change(num, poly_ring_tag)
        den_tag = parent_ring_change(den, poly_ring_tag)
        Q = lcm(Q, den_tag)
        tagged_poly_mqs = num_tag - tag_vars[i] * den_tag
        tagged_mqs[i] = tagged_poly_mqs
    end
    tagged_mqs[end] = Q * sat_var - 1
    # Compute the basis of the MQS in K[T][x][t] such that T < x < t.
    #
    # NOTE: we compute the basis in K[T][x][t], not in K(T)[x][t].
    # This way, we obtain two pieces of information at once:
    # - the kernel of the map from T to x (the algbraic relations of the tags).
    # - the GB of the MQS in K(T)[x][t].
    # ord = Lex()
    ord = DegRevLex([sat_var]) * DegRevLex(orig_vars) * DegRevLex(tag_vars)
    @info """
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
    @info """
    Tagged MQS GB:
    $tagged_mqs_gb
    Relations between tags:
    $relations_between_tags
    """
    # Reduce the fractions with respect to the MQS ideal.
    #
    # NOTE: reduction actually happens in K(T)[x]. So we map polynomials to the
    # parametric ring K(T)[x].
    ring_of_tags, tags = PolynomialRing(K, tag_strings)
    tag_to_gen = Dict(tags[i] => fracs_gen[i] for i in 1:length(fracs_gen))
    if !isempty(intersect(tag_strings, orig_strings))
        @warn """
        There is an intersection between the names of the tag variables and the original variables.
        Tags: $tag_strings
        Original vars: $orig_strings"""
    end
    parametric_ring, _ =
        PolynomialRing(FractionField(ring_of_tags), orig_strings, ordering = :degrevlex)
    relations_between_tags =
        map(poly -> parent_ring_change(poly, ring_of_tags), relations_between_tags)
    param_var_mapping = merge(
        Dict(gens(poly_ring_tag)[2:(nvars(poly_ring) + 1)] .=> gens(parametric_ring)),
        Dict(gens(poly_ring_tag)[(nvars(poly_ring) + 2):end] .=> gens(ring_of_tags)),
    )
    @debug """
    Variable mapping:
    $param_var_mapping
    Parametric ring:
    $parametric_ring
    """
    tagged_mqs_gb_param = map(
        poly -> crude_parent_ring_change(poly, parametric_ring, param_var_mapping),
        tagged_mqs_gb,
    )
    tagged_mqs_gb_param = map(f -> divexact(f, leading_coefficient(f)), tagged_mqs_gb_param)
    @debug "Tagged parametric mqs: $tagged_mqs_gb_param"
    # Reduce each fraction
    var_mapping = Dict(gens(poly_ring) .=> gens(parametric_ring))
    memberships = Vector{Bool}(undef, length(to_be_reduced))
    remainders = Vector{Generic.Frac{T}}(undef, length(to_be_reduced))
    for i in 1:length(to_be_reduced)
        frac = to_be_reduced[i]
        num = crude_parent_ring_change(numerator(frac), parametric_ring, var_mapping)
        den = crude_parent_ring_change(denominator(frac), parametric_ring, var_mapping)
        membership, remainder = check_constructive_field_membership(
            tagged_mqs_gb_param,
            relations_between_tags,
            num,
            den,
        )
        memberships[i] = membership
        remainders[i] = remainder
    end
    return memberships, remainders, relations_between_tags, tag_to_gen
end

"""
    vector_field_along(derivation, directions)

Returns the vector field obtained by applying `derivation` to each element of
`directions`.
"""
function vector_field_along(derivation::Dict{T, U}, directions::AbstractVector) where {T, U}
    new_vector_field =
        Dict{AbstractAlgebra.Generic.Frac{T}, AbstractAlgebra.Generic.Frac{T}}()
    for func in directions
        df = diff_frac(func, derivation)
        new_vector_field[func] = df
    end
    return new_vector_field
end

"""
    reparametrize_with_respect_to(ode, new_states, new_params)

Reparametrizes the `ode` using the given fractional states and parameters.

## Input

- `ode`: an ODE model.
- `new_states`: a vector of new states as fractions in `parent(ode)`.
- `new_params`: a vector of new parameters as fractions in `parent(ode)`.
"""
function reparametrize_with_respect_to(ode, new_states, new_params)
    @assert length(new_states) > 0
    poly_ring = base_ring(parent(first(new_states)))
    # Compute the new dynamics in terms of the original variables.
    # Paying attenton to the order..
    new_vector_field = vector_field_along(ode.x_equations, new_states)
    @info "New vector field:\n$new_vector_field"
    new_states = collect(keys(new_vector_field))
    new_dynamics = [new_vector_field[new_state] for new_state in new_states]
    # Express the new dynamics in terms of new states and new parameters.
    outputs = [ode.y_equations[output] for output in ode.y_vars]
    generating_funcs = vcat(
        new_states,
        new_params,
        ode.u_vars .// one(poly_ring),
        ode.y_vars .// one(poly_ring),
    )
    to_be_reduced_funcs = vcat(new_dynamics, outputs .// one(poly_ring))
    n_active_generators =
        (length(generating_funcs) - length(ode.u_vars) - length(ode.y_vars))
    tag_names = vcat(
        gen_tag_names(n_active_generators, "Internal"),
        gen_tag_names(length(ode.u_vars), "Input"),
        gen_tag_names(length(ode.y_vars), "Output"),
    )
    @info """
    Tag names: 
    $tag_names
    Generating functions:
    $generating_funcs
    To be reduced functions:
    $to_be_reduced_funcs
    """
    membership, new_dynamics_all, implicit_relations, new_vars =
        check_constructive_field_membership(
            RationalFunctionField(generating_funcs),
            to_be_reduced_funcs;
            tag_names = tag_names,
        )
    @assert all(membership)
    ring_of_tags = parent(first(keys(new_vars)))
    tags = gens(ring_of_tags)
    tag_inputs = tags[(n_active_generators + 1):(end - length(ode.y_vars))]
    tag_outputs = tags[(end - length(ode.y_vars) + 1):end]
    new_dynamics_states = new_dynamics_all[1:length(new_states)]
    new_dynamics_outputs = new_dynamics_all[(length(new_states) + 1):end]
    new_outputs = Dict(
        output => dynamic for (output, dynamic) in zip(tag_outputs, new_dynamics_outputs)
    )
    # Old inputs map one to one to new inputs.
    new_inputs = empty(tags)
    if !isempty(ode.u_vars)
        new_inputs = tag_inputs
    end
    @info """
    New state dynamics:
    $new_dynamics_states
    New output dynamics:
    $new_dynamics_outputs
    New inputs:
    $new_inputs"""
    # Construct the new vector field.
    new_vars_vector_field = empty(ode.x_equations)
    for i in 1:length(new_states)
        state = tags[i]
        new_vars_vector_field[state] = new_dynamics_states[i]
    end
    @assert parent(first(keys(new_vars_vector_field))) ==
            base_ring(parent(first(values(new_vars_vector_field)))) ==
            parent(first(keys(new_outputs))) ==
            base_ring(parent(first(values(new_outputs)))) ==
            parent(first(keys(new_vars)))
    @assert base_ring(parent(first(values(new_vars)))) == parent(ode)
    new_vars_vector_field, new_inputs, new_outputs, new_vars, implicit_relations
end

"""
    reparametrize_global(ode, options...)

Casts an incantation and returns a rabbit.

## Options

The function accepts the following optional arguments.

- `seed`: A float in the range from 0 to 1, random seed (default is `seed = 42`). 
- `p`: The probability of correctness (default is `p = 0.99`).
"""
function reparametrize_global(ode::ODE{P}; p = 0.99, seed = 42) where {P}
    Random.seed!(seed)
    id_funcs =
        find_identifiable_functions(ode, with_states = true, simplify = :strong, p = p)
    ode_ring = parent(ode)
    @assert base_ring(parent(first(id_funcs))) == ode_ring
    @info "Constructing a new parametrization"
    contains_states(poly::MPolyElem) = any(x -> degree(poly, x) > 0, ode.x_vars)
    contains_states(func) =
        contains_states(numerator(func)) || contains_states(denominator(func))
    id_funcs_contains_states = filter(contains_states, id_funcs)
    @info """
    Original states: $(ode.x_vars)
    Original params: $(ode.parameters)
    Identifiable and contain states: $(id_funcs_contains_states)"""
    new_states = id_funcs_contains_states
    new_params = setdiff(id_funcs, id_funcs_contains_states)
    @info """
    Reparametrizing with respect to:
    New states: $new_states
    New params: $new_params"""
    new_vector_field, new_inputs, new_outputs, new_vars, implicit_relations =
        reparametrize_with_respect_to(ode, new_states, new_params)
    new_ring = parent(first(keys(new_vector_field)))
    new_ode = ODE{P}(new_vector_field, new_outputs, new_inputs)
    return (new_ode = new_ode, new_vars = new_vars, implicit_relations = implicit_relations)
end
