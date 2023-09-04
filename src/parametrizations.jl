# Reparametrize!

to_fractions(polys::Vector{T}) where {T} = polys .// one(first(polys))
to_fractions(dennums::Vector{Vector{T}}) where {T} =
    StructuralIdentifiability.dennums_to_fractions(dennums)
to_fractions(fracs::Vector{AbstractAlgebra.Generic.Frac{T}}) where {T} = fracs

"""
    check_constructive_field_membership(generators, to_be_reduced)

Returns the unique expression of `to_be_reduced` in terms of the elements of
`generators`.

Follows the vein of Algorithm 1.17 from https://doi.org/10.1006/jsco.1998.0246
"""
function check_constructive_field_membership(generators::AbstractVector, to_be_reduced)
    fracs_gen = to_fractions(generators)
    frac_to_be_reduced = first(to_fractions([to_be_reduced]))
    @assert parent(first(fracs_gen)) == parent(frac_to_be_reduced)
    ring = base_ring(parent(frac_to_be_reduced))
    K = base_ring(ring)
    tag_strings = map(i -> "T$i", 1:length(fracs_gen))
    sat_string = "t"
    @info """
    Tags:
    $(join(map(x -> string(x[1]) * " -> " * string(x[2]),  zip(fracs_gen, tag_strings)), "\t\n"))
    """
    var_strings = vcat(sat_string, map(string, gens(ring)), tag_strings)
    ring_tag, xs_tag = PolynomialRing(K, var_strings, ordering = Nemo.ordering(ring))
    orig_vars = xs_tag[2:(nvars(ring) + 1)]
    tag_vars = xs_tag[(nvars(ring) + 2):end]
    sat_var = xs_tag[1]
    @assert all(<(sat_var), tag_vars)
    @assert all(orig_var -> all(<(orig_var), tag_vars), orig_vars)
    tag_to_gen = Dict(tag_vars[i] => fracs_gen[i] for i in 1:length(fracs_gen))
    @info """
    Original poly ring: $ring
    Tagged poly ring: $ring_tag"""
    tagged_mqs = Vector{elem_type(ring)}()
    num, den = StructuralIdentifiability.unpack_fraction(frac_to_be_reduced)
    # Fraction to be reduced is Num // Den, A is a formal parameter.
    #
    # Construct Num - A * Den to later compute the normal form of it w.r.t. the
    # generators. Then, A = Num / Den would be the desired expression.
    #
    # Note that since A is not present in the generators, the normal form is
    # multiplicative by A, and, thus, can be computed separately for Num and Den.
    to_be_reduced_tag = (
        StructuralIdentifiability.parent_ring_change(num, ring_tag),
        StructuralIdentifiability.parent_ring_change(den, ring_tag),
    )
    Q = one(ring_tag)
    for i in 1:length(fracs_gen)
        num, den = StructuralIdentifiability.unpack_fraction(fracs_gen[i])
        num_tag = StructuralIdentifiability.parent_ring_change(num, ring_tag)
        den_tag = StructuralIdentifiability.parent_ring_change(den, ring_tag)
        Q = lcm(Q, den_tag)
        tagged_poly_mqs = num_tag - tag_vars[i] * den_tag
        push!(tagged_mqs, tagged_poly_mqs)
    end
    push!(tagged_mqs, Q * sat_var - 1)
    # ord = DegRevLex([sat_var]) * DegRevLex(orig_vars) * DegRevLex(tag_vars)
    ord = Lex()
    @info """
    Tagged MQS ideal:
    $tagged_mqs
    Monom ordering:
    $(ord)"""
    tagged_mqs_gb = groebner(tagged_mqs, ordering = ord)
    tags_syzygies = filter(
        poly -> isempty(intersect(vars(poly), vcat(sat_var, orig_vars))),
        tagged_mqs_gb,
    )
    tagged_mqs_gb = setdiff(tagged_mqs_gb, tags_syzygies)
    tagged_mqs_gb = filter(poly -> isempty(intersect(vars(poly), [sat_var])), tagged_mqs_gb)
    @info """
    Tagged MQS GB:
    $tagged_mqs_gb
    Syzygies of tags:
    $tags_syzygies
    To be reduced:
    $to_be_reduced_tag
    """
    function switch_ring_to_parametric(poly, new_ring)
        param_ring = base_ring(base_ring(new_ring))
        params = gens(param_ring)
        orig_vars = gens(ring)
        new_poly = zero(new_ring)
        for (i, term) in enumerate(terms(poly))
            new_coeff = one(param_ring) * coeff(poly, i)
            new_monom = one(new_ring)
            for var in vars(term)
                exp = degree(term, var)
                if string(var) in map(string, params)
                    new_coeff *=
                        StructuralIdentifiability.parent_ring_change(var, param_ring)^exp
                else
                    new_monom *=
                        StructuralIdentifiability.parent_ring_change(var, new_ring)^exp
                end
            end
            new_poly += new_coeff * new_monom
        end
        new_poly
    end
    ring_of_tags, = PolynomialRing(K, tag_strings)
    parametric_ring, _ =
        PolynomialRing(FractionField(ring_of_tags), map(string, orig_vars), ordering = :lex)
    tagged_mqs_gb =
        map(poly -> switch_ring_to_parametric(poly, parametric_ring), tagged_mqs_gb)
    tags_syzygies =
        map(poly -> switch_ring_to_parametric(poly, parametric_ring), tagged_mqs_gb)
    to_be_reduced_tag =
        map(poly -> switch_ring_to_parametric(poly, parametric_ring), to_be_reduced_tag)
    _, num_rem = divrem(to_be_reduced_tag[1], tagged_mqs_gb)
    _, den_rem = divrem(to_be_reduced_tag[2], tagged_mqs_gb)
    @info "" num_rem den_rem
    remainder = num_rem // den_rem
    _, num_factored = divrem(numerator(remainder), tags_syzygies)
    _, den_factored = divrem(denominator(remainder), tags_syzygies)
    if iszero(den_factored) ||
       !isempty(
           intersect(
               map(string, vars(num_factored)),
               map(string, vcat(sat_var, orig_vars)),
           ),
       ) ||
       !isempty(
           intersect(
               map(string, vars(den_factored)),
               map(string, vcat(sat_var, orig_vars)),
           ),
       )
        @warn """
        The fraction ($(frac_to_be_reduced)) is not a function of the generators ($(fracs_gen)).
        Normal form: $remainder
        Normal form num (syzygies factored out): $num_factored
        Normal form den (syzygies factored out): $den_factored
        """
    end
    num_factored =
        StructuralIdentifiability.parent_ring_change(coeff(num_factored, 1), ring_tag)
    den_factored =
        StructuralIdentifiability.parent_ring_change(coeff(den_factored, 1), ring_tag)
    remainder = num_factored // den_factored
    remainder, tag_to_gen
end

# Same as above, but reduces multiple fractions at once.
function check_constructive_field_membership(
    generators::AbstractVector,
    to_be_reduced::AbstractVector,
)
    fracs_gen = to_fractions(generators)
    fracs_to_be_reduced = to_fractions(to_be_reduced)
    T = elem_type(base_ring(parent(first(fracs_to_be_reduced))))
    remainders = Vector{Generic.Frac{T}}(undef, length(fracs_to_be_reduced))
    tag_to_gen = Dict{T, Generic.Frac{T}}()
    for i in 1:length(fracs_to_be_reduced)
        frac = fracs_to_be_reduced[i]
        remainder, tag_to_gen = check_constructive_field_membership(fracs_gen, frac)
        remainders[i] = remainder
    end
    @assert length(unique(parent, remainders)) == 1
    remainders, tag_to_gen
end

"""
    vector_field_along(derivation, directions)

Returns the vector field obtained by applying `derivation` to each element of
`directions`.
"""
function vector_field_along(derivation::Dict{T, U}, directions::AbstractVector) where {T, U}
    fracs = to_fractions(directions)
    new_vector_field =
        Dict{AbstractAlgebra.Generic.Frac{T}, AbstractAlgebra.Generic.Frac{T}}()
    for func in fracs
        df = StructuralIdentifiability.diff_frac(func, derivation)
        new_vector_field[func] = df
    end
    new_vector_field
end

"""
    reparametrize_with_respect_to(ode, new_states, new_params)

Reparametrizes the `ode` using the given states and parameters.

## Input

- `ode`: an ODE model.
- `new_states`: a vector of new states as functions in `parent(ode)`.
- `new_params`: a vector of new parameters as functions in `parent(ode)`.
"""
function reparametrize_with_respect_to(ode, new_states, new_params)
    @assert length(new_states) + length(new_params) > 0
    # Compute the new dynamics in terms of the original variables.
    # Paying attenton to the order..
    new_vector_field = vector_field_along(ode.x_equations, new_states)
    states = collect(keys(new_vector_field))
    dynamics = [new_vector_field[state] for state in states]
    # Express the new dynamics in terms of new states and new parameters.
    generating_funcs = vcat(states, new_params, ode.u_vars)
    new_vars_dynamics, new_vars =
        check_constructive_field_membership(generating_funcs, dynamics)
    tag_ring = parent(first(keys(new_vars)))
    # Express the existing outputs in terms of new states and new parameters.
    outputs = ode.y_vars
    new_outputs_dynamics, _ = check_constructive_field_membership(
        generating_funcs,
        [ode.y_equations[output] for output in outputs],
    )
    new_outputs = Dict(
        StructuralIdentifiability.parent_ring_change(output, tag_ring) => dynamic for
        (output, dynamic) in zip(outputs, new_outputs_dynamics)
    )
    # NOTE: old inputs map one to one to new inputs.
    inputs = ode.u_vars
    if !isempty(inputs)
        new_inputs_dynamics, _ =
            check_constructive_field_membership(generating_funcs, inputs)
        new_inputs = Dict(
            StructuralIdentifiability.parent_ring_change(input, tag_ring) => new_input
            for (input, new_input) in zip(inputs, new_inputs_dynamics)
        )
    else
        new_inputs = empty(new_outputs)
    end
    # Construct the new vector field.
    new_vars_vector_field = empty(ode.x_equations)
    state_to_new_var = Dict(v => k for (k, v) in new_vars)
    for i in 1:length(states)
        state = states[i]
        new_vars_vector_field[state_to_new_var[state]] = new_vars_dynamics[i]
    end
    @assert parent(first(keys(new_vars_vector_field))) ==
            base_ring(parent(first(values(new_vars_vector_field)))) ==
            parent(first(keys(new_outputs))) ==
            base_ring(parent(first(values(new_outputs)))) ==
            parent(first(keys(new_vars)))
    @assert base_ring(parent(first(values(new_vars)))) == parent(ode)
    new_vars_vector_field, new_inputs, new_outputs, new_vars
end

function reparametrize_global(ode::StructuralIdentifiability.ODE{P}) where {P}
    id_funcs = StructuralIdentifiability.find_identifiable_functions(
        ode,
        with_states = true,
        strategy = (:hybrid,),
    )
    @assert base_ring(parent(first(id_funcs))) == parent(ode)
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
    new_vector_field, new_inputs, new_outputs, new_vars =
        reparametrize_with_respect_to(ode, new_states, new_params)
    new_ring = parent(first(keys(new_vector_field)))
    new_vars_trimmed = union(map(vars, collect(keys(new_vector_field)))...)
    new_vars_trimmed =
        union(new_vars_trimmed, map(vars, collect(values(new_vector_field)))...)
    # new_vars_trimmed = union(new_vars_trimmed, map(vars, collect(keys(new_inputs)))...)
    new_vars_trimmed = union(new_vars_trimmed, map(vars, collect(keys(new_outputs)))...)
    new_vars_trimmed = union(new_vars_trimmed, map(vars, collect(values(new_outputs)))...)
    new_ring_trimmed, new_vars_trimmed = PolynomialRing(
        base_ring(new_ring),
        map(string, new_vars_trimmed),
        ordering = Nemo.ordering(new_ring),
    )
    new_ode = StructuralIdentifiability.ODE{P}(
        Dict(
            StructuralIdentifiability.parent_ring_change(k, new_ring_trimmed) =>
                StructuralIdentifiability.parent_ring_change(v, new_ring_trimmed) for
            (k, v) in new_vector_field
        ),
        Dict(
            StructuralIdentifiability.parent_ring_change(k, new_ring_trimmed) =>
                StructuralIdentifiability.parent_ring_change(v, new_ring_trimmed) for
            (k, v) in new_outputs
        ),
        map(
            f -> StructuralIdentifiability.parent_ring_change(
                numerator(f),
                new_ring_trimmed,
            ),
            collect(values(new_inputs)),
        ),
    )
    new_vars = Dict(
        StructuralIdentifiability.parent_ring_change(k, new_ring_trimmed) => v for
        (k, v) in new_vars
    )
    @assert base_ring(parent(first(values(new_vars)))) == parent(ode)
    return new_ode, new_vars
end
