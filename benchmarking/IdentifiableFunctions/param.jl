# This is a proof of concept implementation of global model reparametrization.

using StructuralIdentifiability
using AbstractAlgebra, Nemo, Groebner

to_fractions(polys::Vector{T}) where {T} = polys .// one(first(polys))
to_fractions(dennums::Vector{Vector{T}}) where {T} =
    StructuralIdentifiability.dennums_to_fractions(dennums)
to_fractions(fracs::Vector{AbstractAlgebra.Generic.Frac{T}}) where {T} = fracs

"""
    check_constructive_field_membership(generators, to_be_reduced)

Returns the unique expression of `to_be_reduced` in terms of the elements of
`generators`.

Assumes that `generators` are algebraically independent.

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
    $(join(map(x -> string(x[1]) * " -> " * string(x[2]), zip(fracs_gen, tag_strings)), "\t\n"))
    """
    var_strings = vcat(sat_string, map(string, gens(ring)), tag_strings)
    ring_tag, xs_tag = polynomial_ring(K, var_strings, ordering = Nemo.ordering(ring))
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
    # ord = DegRevLex([sat_var]) * DegRevLex(orig_vars) * Groebner.DegRevLex(tag_vars)
    ord = Lex()
    tagged_mqs_gb = groebner(tagged_mqs, ordering = ord)
    @info """
    Tagged MQS ideal:
    $tagged_mqs
    Tagged MQS GB:
    $tagged_mqs_gb
    GB ordering:
    $(ord)
    To be reduced:
    $to_be_reduced_tag
    """
    _, num_rem = divrem(to_be_reduced_tag[1], tagged_mqs_gb)
    _, den_rem = divrem(to_be_reduced_tag[2], tagged_mqs_gb)
    remainder = num_rem // den_rem
    if !all(orig_var -> degree(numerator(remainder), orig_var) == 0, orig_vars) ||
            !all(orig_var -> degree(denominator(remainder), orig_var) == 0, orig_vars)
        @warn """
        The fraction ($(frac_to_be_reduced)) is not a function of the generators ($(fracs_gen)).
        Normal form: $remainder
        """
    end
    return remainder, tag_to_gen
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
    return remainders, tag_to_gen
end

#####################
# Sanity test!
R, (a, b) = Nemo.QQ["a", "b"]

fracs_generators = [a^2, (a + b) // (b)]
to_be_reduced = [a^4, (a + b - 5a^2 * b) // (b), (a + b) // (b * a^20)]

rem_tags, tag_to_gen = check_constructive_field_membership(fracs_generators, to_be_reduced)
@info "" rem_tags tag_to_gen
#=
┌ Info: 
│   rem_tags =
│    3-element Vector{AbstractAlgebra.Generic.Frac{QQMPolyRingElem}}:
│     T1^2
│     -5*T1 + T2
│     T2//T1^10
│   tag_to_gen =
│    Dict{QQMPolyRingElem, AbstractAlgebra.Generic.Frac{QQMPolyRingElem}} with 2 entries:
│      T1 => a^2
└      T2 => (a + b)//b
=#
#####################

"""
    vector_field_along(derivation, directions)

Returns the vector field obtained by applying `derivation` to each element of
`directions`.
"""
function vector_field_along(derivation, directions::AbstractVector)
    fracs = to_fractions(directions)
    new_vector_field = Dict()
    for func in fracs
        df = StructuralIdentifiability.diff_frac(func, derivation)
        new_vector_field[func] = df
    end
    return new_vector_field
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
    # Compute the new dynamics in terms of the original variables.
    # Paying attenton to the order..
    new_vector_field = vector_field_along(ode.x_equations, new_states)
    states = collect(keys(new_vector_field))
    dynamics = [new_vector_field[state] for state in states]
    # Express the new dynamics in terms of new states and new parameters.
    generating_funcs = vcat(states, new_params)
    new_vars_dynamics, new_vars =
        check_constructive_field_membership(generating_funcs, dynamics)
    # Express the existing outputs in terms of new states and new parameters.
    outputs = ode.y_vars
    new_outputs_dynamics, _ = check_constructive_field_membership(
        generating_funcs,
        [ode.y_equations[output] for output in outputs],
    )
    new_outputs =
        Dict(output => dynamic for (output, dynamic) in zip(outputs, new_outputs_dynamics))
    # Construct the new vector field.
    new_vars_vector_field = empty(new_vector_field)
    state_to_new_var = Dict(v => k for (k, v) in new_vars)
    for i in 1:length(states)
        state = states[i]
        new_vars_vector_field[state_to_new_var[state]] = new_vars_dynamics[i]
    end
    return new_vars_vector_field, new_outputs, new_vars
end

#####################
# Sanity test!
# 1.
ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = x1 + x2 + a + b,
    x2'(t) = x1 + x2,
    y(t) = x1 + x2
)

new_vector_field, new_outputs, new_vars =
    reparametrize_with_respect_to(ode, [x1 + x2], [a + b])
@info "" new_vector_field new_outputs new_vars
#=
┌ Info: 
│   new_vector_field =
│    Dict{Any, Any} with 1 entry:
│      T1 => 2*T1 + T2
│   new_outputs =
│    Dict{QQMPolyRingElem, AbstractAlgebra.Generic.Frac{QQMPolyRingElem}} with 1 entry:
│      y => T1
│   new_vars =
│    Dict{QQMPolyRingElem, AbstractAlgebra.Generic.Frac{QQMPolyRingElem}} with 2 entries:
│      T2 => a + b
└      T1 => x2 + x1
=#

# 2.
