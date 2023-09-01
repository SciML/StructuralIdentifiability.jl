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

R, (x1, x2, a, b) = Nemo.QQ["x1", "x2", "a", "b"]

f = [a // b, b^2]
rem_tag, tag_to_gen = tagged_normal_form(f, a^2)
