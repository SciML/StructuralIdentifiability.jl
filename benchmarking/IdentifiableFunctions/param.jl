using AbstractAlgebra, Nemo, Groebner

R, (x1, x2, a, b) = Nemo.QQ["x1", "x2", "a", "b"]

f = [x1 + x2, a + b]
rff = StructuralIdentifiability.RationalFunctionField(f)

StructuralIdentifiability.simplified_generating_set(rff)

to_fractions(polys::Vector{T}) where {T} = polys .// one(first(polys))
to_fractions(dennums::Vector{Vector{T}}) where {T} =
    StructuralIdentifiability.dennums_to_fractions(dennums)
to_fractions(fracs::Vector{AbstractAlgebra.Generic.Frac{T}}) where {T} = fracs

function tagged_normal_form(generators, to_be_reduced)
    fracs = to_fractions(generators)
    to_be_reduced_frac = first(to_fractions([to_be_reduced]))
    @assert parent(first(fracs)) == parent(to_be_reduced_frac)
    ring = base_ring(parent(to_be_reduced_frac))
    K = base_ring(ring)
    tag_strings = map(i -> "T$i", 1:length(fracs))
    sat_string = "t"
    @info """
    Tags:
    $(join(map(x -> string(x[1]) * " -> " * string(x[2]),  zip(fracs, tag_strings)), "\t\n"))
    """
    var_strings = vcat(map(string, gens(ring)), tag_strings, sat_string)
    ring_tag, xs_tag = PolynomialRing(K, var_strings, ordering = Nemo.ordering(ring))
    tag_vars = xs_tag[(nvars(ring) + 1):(end - 1)]
    sat_var = xs_tag[end]
    @info """
    Original poly ring: $ring
    Tagged poly ring: $ring_tag"""
    tagged_mqs = Vector{elem_type(ring)}()
    num, den = StructuralIdentifiability.unpack_fraction(to_be_reduced_frac)
    # Num - A * Den
    to_be_reduced_tag = (
        StructuralIdentifiability.parent_ring_change(num, ring_tag, matching = :byname),
        StructuralIdentifiability.parent_ring_change(den, ring_tag, matching = :byname),
    )
    Q = one(ring_tag)
    for i in 1:length(fracs)
        num, den = StructuralIdentifiability.unpack_fraction(fracs[i])
        num_tag =
            StructuralIdentifiability.parent_ring_change(num, ring_tag, matching = :byname)
        den_tag =
            StructuralIdentifiability.parent_ring_change(den, ring_tag, matching = :byname)
        Q = lcm(Q, den_tag)
        tagged_poly_mqs = num_tag - tag_vars[i] * den_tag
        push!(tagged_mqs, tagged_poly_mqs)
    end
    push!(tagged_mqs, Q * sat_var - 1)
    tagged_mqs_gb = Groebner.groebner(tagged_mqs)
    @info """
    Tagged MQS ideal:
    $tagged_mqs
    Tagged MQS GB:
    $tagged_mqs_gb
    GB ordering:
    $(Nemo.ordering(parent(first(tagged_mqs_gb))))
    To be reduced:
    $to_be_reduced_tag
    """
    _, num_rem = divrem(to_be_reduced_tag[1], tagged_mqs_gb)
    _, den_rem = divrem(to_be_reduced_tag[2], tagged_mqs_gb)
    # A = Num / Den
    remainder = num_rem // den_rem
    tag_to_gen = Dict(tag_vars[i] => fracs[i] for i in 1:length(fracs))
    remainder, tag_to_gen
end

function vector_field_along(derivation, directions)
    fracs = to_fractions(directions)
    new_vector_field = Dict()
    for func in fracs
        df = StructuralIdentifiability.diff_frac(func, derivation)
        new_vector_field[func] = df
    end
    new_vector_field
end

function reparametrize(ode)
    id_funcs = StructuralIdentifiability.find_identifiable_functions(
        ode,
        with_states = true,
        strategy = (:normalforms, 3),
    )
    @info "Constructing a new parametrization"
    poly_ring = base_ring(parent(id_funcs[1]))
    # poly_ring, _ = PolynomialRing(
    #     base_ring(poly_ring),
    #     vcat(map(string, gens(poly_ring)), map(string, ode.u_vars)),
    #     ordering = Nemo.ordering(poly_ring),
    # )
    reparametrize_with_respect_to(ode, new_states, new_params)
    states_id_ring =
        map(x -> StructuralIdentifiability.parent_ring_change(x, poly_ring), ode.x_vars)
    derivation_id_ring = Dict(
        StructuralIdentifiability.parent_ring_change(k, poly_ring) =>
            StructuralIdentifiability.parent_ring_change(v, poly_ring) for
        (k, v) in ode.x_equations
    )
    contains_states(poly::MPolyElem) = any(x -> degree(poly, x) > 0, states_id_ring)
    contains_states(func) =
        contains_states(numerator(func)) || contains_states(denominator(func))
    id_funcs_contains_states = filter(contains_states, id_funcs)
    @info """
    Original states: $(ode.x_vars)
    Original params: $(ode.parameters)
    Identifiable states: $(id_funcs_contains_states)
    Derivation: $(derivation_id_ring)"""
    new_vector_field = vector_field_along(derivation_id_ring, id_funcs_contains_states)
    T = elem_type(poly_ring)
    new_tagged_vector_field = Dict{T, AbstractAlgebra.Generic.Frac{T}}()
    new_vars = Dict{T, AbstractAlgebra.Generic.Frac{T}}()
    id_func_to_index = Dict(id_funcs[i] => i for i in 1:length(id_funcs))
    for (func, df) in new_vector_field
        rem_tag, tag_to_gen = tagged_normal_form(id_funcs, df)
        tags = gens(base_ring(parent(rem_tag)))
        new_tagged_vector_field[tags[id_func_to_index[func] + nvars(poly_ring)]] = rem_tag
        for (k, v) in tag_to_gen
            new_vars[k] = v
        end
    end
    new_poly_ring, new_xs = PolynomialRing(
        base_ring(poly_ring),
        vcat(
            map(string, gens(parent(first(keys(new_vars)))))[(nvars(poly_ring) + 1):(end - 1)],
            map(string, ode.y_vars),
            # map(string, ode.u_vars),
        ),
        ordering = Nemo.ordering(poly_ring),
    )
    new_ys = new_xs[(length(new_xs) - length(ode.y_vars) + 1):(end)]
    new_outputs = Dict{T, AbstractAlgebra.Generic.Frac{T}}()
    for (i, (k, v)) in enumerate(ode.y_equations)
        y_dynamics = StructuralIdentifiability.parent_ring_change(v, poly_ring)
        rem_tag, tag_to_gen = tagged_normal_form(id_funcs, y_dynamics)
        new_outputs[new_ys[i]] =
            StructuralIdentifiability.parent_ring_change(rem_tag, new_poly_ring)
    end
    new_vector_field = Dict(
        StructuralIdentifiability.parent_ring_change(k, new_poly_ring) =>
            StructuralIdentifiability.parent_ring_change(v, new_poly_ring) for
        (k, v) in new_tagged_vector_field
    )
    new_vars = Dict(
        StructuralIdentifiability.parent_ring_change(k, new_poly_ring) =>
            StructuralIdentifiability.parent_ring_change(v, parent(ode)) for
        (k, v) in new_vars
    )
    for (i, new_y) in enumerate(new_ys)
        new_vars[new_y] = ode.y_vars[i] // one(parent(ode))
    end
    @assert parent(first(keys(new_vector_field))) ==
            base_ring(parent(first(values(new_vector_field)))) ==
            parent(first(keys(new_outputs))) ==
            base_ring(parent(first(values(new_outputs)))) ==
            parent(first(keys(new_vars)))
    @assert base_ring(parent(first(values(new_vars)))) == parent(ode)
    new_vector_field, new_outputs, new_vars
end

ode = StructuralIdentifiability.@ODEmodel(
    x1'(t) = x1 + x2 + a + b,
    x2'(t) = x1 + x2 + u,
    y(t) = x1 + x2
)

StructuralIdentifiability.find_identifiable_functions(ode, with_states = true)

new_vector_field, new_outputs, new_vars = reparametrize(ode)

new_ode = StructuralIdentifiability.ODE{keytype(new_vector_field)}(
    new_vector_field,
    new_outputs,
    Vector{P}(),
)

StructuralIdentifiability.assess_identifiability(new_ode)

R, (x1, x2, a, b) = Nemo.QQ["x1", "x2", "a", "b"]

f = [a]
rem_tag, tag_to_gen = tagged_normal_form(f, b)

rem_tag, tag_to_gen = tagged_normal_form(f, b)

hiv = StructuralIdentifiability.@ODEmodel(
    w'(t) = -b * w(t) + c * w(t) * x(t) * y(t) - c * w(t) * q * y(t),
    v'(t) = k * y(t) - v(t) * u,
    x'(t) = lm - x(t) * d - x(t) * v(t) * beta,
    z'(t) = c * w(t) * q * y(t) - h * z(t),
    y'(t) = x(t) * v(t) * beta - a * y(t),
    y2(t) = z(t),
    y1(t) = w(t)
)

StructuralIdentifiability.assess_identifiability(hiv)

new_vector_field, new_outputs, new_vars = reparametrize(hiv)

sirs = StructuralIdentifiability.@ODEmodel(
    s'(t) = mu - mu * s(t) - b0 * (1 + b1 * x1(t)) * i(t) * s(t) + g * r(t),
    i'(t) = b0 * (1 + b1 * x1(t)) * i(t) * s(t) - (nu + mu) * i(t),
    r'(t) = nu * i(t) - (mu + g) * r(t),
    x1'(t) = -M * x2(t),
    x2'(t) = M * x1(t),
    y1(t) = i(t),
    y2(t) = r(t)
)
vector_field_along(sirs.x_equations, [x1 * x2])

new_vector_field, new_outputs, new_vars = reparametrize(sirs)

new_ode = StructuralIdentifiability.ODE{keytype(new_vector_field)}(
    new_vector_field,
    new_outputs,
    Vector{P}(),
)

StructuralIdentifiability.assess_identifiability(sirs)
StructuralIdentifiability.assess_identifiability(new_ode)

funcs = StructuralIdentifiability.find_identifiable_functions(
    sirs,
    with_states = true,
    strategy = (:normalforms, 5),
)
