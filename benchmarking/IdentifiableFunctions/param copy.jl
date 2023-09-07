# This is a proof of concept implementation of global model reparametrization.

# using StructuralIdentifiability
using AbstractAlgebra, Nemo, Groebner

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

R, (a, b) = Nemo.QQ["a", "b"]

T1, T2, T3 = a^2 + b^2, a^3 + b^3, a^4 + b^4
T1^6 - 4 * T1^3 * T2^2 - 3 * T1^2 * T3^2 + 12 * T1 * T2^2 * T3 - 4 * T2^4 - 2 * T3^3

fracs_generators = [a^2 + b^2]
to_be_reduced = [a^2 + b^2]
rem_tags, tag_to_gen = check_constructive_field_membership(fracs_generators, to_be_reduced)
@info "" rem_tags tag_to_gen

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
│    3-element Vector{AbstractAlgebra.Generic.Frac{fmpq_mpoly}}:
│     T1^2
│     -5*T1 + T2
│     T2//T1^10
│   tag_to_gen =
│    Dict{fmpq_mpoly, AbstractAlgebra.Generic.Frac{fmpq_mpoly}} with 2 entries:
│      T1 => a^2
└      T2 => (a + b)//b
=#
#####################

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

##########################

covid = StructuralIdentifiability.@ODEmodel(
    S'(t) = -b * S(t) * (I(t) + J(t) + q * A(t)) * Ninv(t),
    E'(t) = b * S(t) * (I(t) + J(t) + q * A(t)) * Ninv(t) - k * E(t),
    A'(t) = k * (1 - r) * E(t) - g1 * A(t),
    I'(t) = k * r * E(t) - (alpha + g1) * I(t),
    J'(t) = alpha * I(t) - g2 * J(t),
    C'(t) = alpha * I(t),
    Ninv'(t) = 0,
    y(t) = C(t),
    y2(t) = Ninv(t)
)

id_funcs = StructuralIdentifiability.find_identifiable_functions(
    covid,
    with_states = true,
    strategy = (:hybrid,),
)

new_states = [J, C, I, r * E, r * S, q * A, A // (r * E - E)]
new_params = [Ninv, g1, k, g2, alpha, b]

new_vector_field, new_inputs, new_outputs, new_vars =
    reparametrize_with_respect_to(covid, new_states, new_params)
@info "" new_vector_field new_inputs new_outputs new_vars

##########################

ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = a * x1 - b * x1 * x2 + u(t),
    x2'(t) = -c * x2 + d * x1 * x2,
    y(t) = x1
)

StructuralIdentifiability.find_identifiable_functions(
    ode,
    with_states = true,
    strategy = (:hybrid,),
)

new_vector_field, new_inputs, new_outputs, new_vars =
    reparametrize_with_respect_to(ode, [x1, b * x2], [a, c, d])
@info "" new_vector_field new_inputs new_outputs new_vars

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
│    Dict{fmpq_mpoly, AbstractAlgebra.Generic.Frac{fmpq_mpoly}} with 1 entry:
│      y => T1
│   new_vars =
│    Dict{fmpq_mpoly, AbstractAlgebra.Generic.Frac{fmpq_mpoly}} with 2 entries:
│      T2 => a + b
└      T1 => x2 + x1
=#

# 2.
#! format: off
ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = a * x1, 
    x2'(t) = b * x2, 
    y(t) = x1*x2
)

id_funcs = StructuralIdentifiability.find_identifiable_functions(
    ode,
    with_states = true,
    strategy = (:hybrid,),
)
@info "" id_funcs
#=
┌ Info: 
│   id_funcs =
│    2-element Vector{AbstractAlgebra.Generic.Frac{fmpq_mpoly}}:
│     x2*x1
└     a + b
=#

new_states = [x1 * x2]
new_params = [a + b]

new_vector_field, new_outputs, new_vars = reparametrize_with_respect_to(ode, new_states, new_params)
@info "" new_vector_field new_outputs new_vars
#=
┌ Info: 
│   new_vector_field =
│    Dict{Any, Any} with 1 entry:
│      T1 => T1*T2
│   new_outputs =
│    Dict{fmpq_mpoly, AbstractAlgebra.Generic.Frac{fmpq_mpoly}} with 1 entry:
│      y => T1
│   new_vars =
│    Dict{fmpq_mpoly, AbstractAlgebra.Generic.Frac{fmpq_mpoly}} with 2 entries:
│      T2 => a + b
└      T1 => x2*x1
=#

# 3.
#! format: off
ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = a * x1, 
    x2'(t) = b * x2, 
    y(t) = x1 + x2
)

id_funcs = StructuralIdentifiability.find_identifiable_functions(
    ode,
    with_states = true,
    strategy = (:hybrid,),
)

T1 = -a*x2 + a*x1 + b*x2 - b*x1
T2 = x2 + x1
T3 = x2*x1
T4 = a*b
T5 = a + b

T1^2 + 4*T2^2*T4 - T2^2*T5^2 - 16*T3*T4 + 4*T3*T5^2

new_vector_field, new_inputs, new_outputs, new_vars = StructuralIdentifiability.reparametrize_with_respect_to(
    ode, [x2 + x1, -a*x2 + a*x1 + b*x2 - b*x1], [a*b, a + b]
)

ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = a * x1 - b * x1 * x2 + u(t),
    x2'(t) = -c * x2 + d * x1 * x2,
    y(t) = x1
)
id_funcs = StructuralIdentifiability.find_identifiable_functions(
    ode,
    with_states = true,
    strategy = (:hybrid,),
)

new_ode, new_vars = StructuralIdentifiability.reparametrize_global(ode)

T3 = x2*x1
T4= a*b
T2= x2 + x1
T5= a + b
T1= -a*x2 + a*x1 + b*x2 - b*x1

id_funcs = StructuralIdentifiability.find_identifiable_functions(
    ode,
    with_states = true,
    strategy = (:hybrid,),
)
@info "" id_funcs
#=
┌ Info: 
│   id_funcs =
│    5-element Vector{AbstractAlgebra.Generic.Frac{fmpq_mpoly}}:
│     x2*x1
│     a*b
│     x2 + x1
│     a + b
└     a*x2 - a*x1 - b*x2 + b*x1
=#

new_states = [x1 * x2, x1 + x2, a*x2 - a*x1 - b*x2 + b*x1]
new_params = [a + b, a*b]

new_vector_field, new_iutputs, new_vars = reparametrize_with_respect_to(ode, new_states, new_params)
@info "" new_vector_field new_outputs new_vars
#=
┌ Info: 
│   new_vector_field =
│    Dict{Any, Any} with 3 entries:
│      T1 => 1//2*T1*T4 - 1//2*T3
│      T2 => T2*T4
│      T3 => -1//2*T1*T4^2 + 2*T1*T5 + 1//2*T3*T4
│   new_outputs =
│    Dict{fmpq_mpoly, AbstractAlgebra.Generic.Frac{fmpq_mpoly}} with 1 entry:
│      y => T1
│   new_vars =
│    Dict{fmpq_mpoly, AbstractAlgebra.Generic.Frac{fmpq_mpoly}} with 5 entries:
│      T1 => x2 + x1
│      T2 => x2*x1
│      T3 => a*x2 - a*x1 - b*x2 + b*x1
│      T4 => a + b
└      T5 => a*b
=#
#! format: on
#####################

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

#####################

using StructuralIdentifiability

# Lotka-Volterra,
# parameter b and state x2 are not identifiabile
ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = a * x1 - b * x1 * x2 + u(t),
    x2'(t) = -c * x2 + d * x1 * x2,
    x3'(t) = x3,
    x4'(t) = x1,
    y(t) = x1
)

new_ode, new_vars = reparametrize_global(ode)
@info "" new_ode new_vars
#=
new_ode =
    T2'(t) = -T2(t)*T1(t) + T2(t)*T5 + T6(t)
    T1'(t) = T2(t)*T1(t)*T4 - T1(t)*T3
    y(t) = T2(t)

new_vars =
Dict{Nemo.fmpq_mpoly, AbstractAlgebra.Generic.Frac{Nemo.fmpq_mpoly}} with 6 entries:
    T2 => x1
    T3 => c
    T1 => b*x2
    T5 => a
    T6 => u
    T4 => d
=#

assess_identifiability(new_ode)
#=
Dict{Any, Symbol} with 5 entries:
  T2 => :globally
  T3 => :globally
  T1 => :globally
  T5 => :globally
  T4 => :globally
=#
