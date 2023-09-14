
# Maintans a row echelon form of a set of vectors over the integrals.
# Works well when the ambient dimension is small.
mutable struct TinyRowEchelonForm{T}
    rows::Vector{Vector{T}}
    pivot_cols::Vector{Int}

    function TinyRowEchelonForm{T}() where {T}
        new(Vector{Vector{T}}(), Vector{Int}())
    end
end

function Base.push!(tref::TinyRowEchelonForm{T}, vect::Vector{T}) where {T}
    @assert count(!iszero, vect) == 1
    @assert !in(tref, vect)
    pivot = findfirst(!iszero, vect)
    push!(tref.rows, vect)
    push!(tref.pivot_cols, pivot)
    return vect
end

function Base.in(tref::TinyRowEchelonForm{T}, vect::Vector{T}) where {T}
    nnz_inds = findall(!iszero, vect)
    for i in nnz_inds
        rref_idx = findfirst(pivot_col -> pivot_col == i, tref.pivot_cols)
        isnothing(rref_idx) && return false
        lead_vect = vect[i]
        lead_rref = tref.rows[rref_idx][i]
        flag, _ = divides(lead_vect, lead_rref)
        !flag && return false
    end
    return true
end

# ------------------------------------------------------------------------------

"""
    local_normal_forms(mqs, field, up_to_degree)

Computes the normal forms of MQS-monomials modulo the MQS-ideal `mqs`
specialized at a random point. Considers monomials up to the total degree
`up_to_degree` over the given `field`.

Ignores any monomials whose exponent vectors are present in `stop_vectors`.

Returns a tuple (`normal_forms`, `monomials`).
"""
function local_normal_forms(
    mqs::IdealMQS,
    finite_field,
    up_to_degree::Integer;
    stop_vectors = TinyRowEchelonForm{Int}(),
)
    ring_param = ParamPunPam.parent_params(mqs)
    point_ff = ParamPunPam.distinct_nonzero_points(finite_field, nvars(ring_param))
    point_ff_ext = vcat(point_ff, one(finite_field))
    gens_ff_spec = specialize_mod_p(mqs, point_ff)
    gb_ff_spec = Groebner.groebner(gens_ff_spec)
    ring_ff = parent(gb_ff_spec[1])
    xs_ff = gens(ring_ff)
    normal_forms_ff = Vector{elem_type(ring_ff)}(undef, 0)
    monoms_ff = Vector{elem_type(ring_ff)}(undef, 0)
    @assert mqs.sat_var_index == length(xs_ff)
    xs_ff = xs_ff[1:(end - 1)]
    pivot_vectors = map(f -> exponent_vector(f, 1), xs_ff)
    @debug """
    variables in the finite field: $(xs_ff)
    gb parent: $(ring_ff)
    specialized gb: $(gb_ff_spec)
    Evaluation point: $point_ff_ext"""
    # Comput normal forms of all possible monomials of degrees from `1` to
    # `up_to_degree`
    for deg in 1:up_to_degree
        for combination in Combinatorics.with_replacement_combinations(pivot_vectors, deg)
            exp_vect = sum(combination)
            if in(stop_vectors, exp_vect)
                @debug "Skipping exponent vector $exp_vect"
                continue
            end
            monom_ff = ring_ff([one(finite_field)], [exp_vect])
            monom_ff_spec = evaluate(monom_ff, point_ff_ext)
            monom_mqs_ff_spec = monom_ff - monom_ff_spec
            divisors_ff, nf_ff = divrem(monom_mqs_ff_spec, gb_ff_spec)
            @debug """
            The normal form of $monom_mqs_ff_spec is:
            normalform = $nf_ff
            divisors = $divisors_ff"""
            push!(monoms_ff, monom_ff)
            push!(normal_forms_ff, nf_ff)
        end
    end
    return normal_forms_ff, monoms_ff
end

# Linearly reduce polys[index] w.r.t polys[1..index-1].
# Return the result and divisors.
function reduce_step_forward(polys, lead_monoms, index, field)
    fi = polys[index]
    multipliers = Vector{Tuple{Int, elem_type(field)}}()
    @inbounds for j in (index - 1):-1:1
        iszero(fi) && break
        fj = polys[j]
        iszero(fj) && continue
        leadj = lead_monoms[j]
        ci = coeff(fi, leadj)
        # If fi contains the lead of fj
        iszero(ci) && continue
        cj = leading_coefficient(fj)
        cij = div(ci, cj)
        fi = fi - cij * fj
        push!(multipliers, (j, -cij))
    end
    return fi, multipliers
end

# Reduce polys[1..index-1] w.r.t polys[index] inplace.
function reduce_step_backward!(polys, lead_monoms, index, field)
    fi = polys[index]
    @assert !iszero(fi)
    leadi = lead_monoms[index]
    ci = leading_coefficient(fi)
    multipliers = Vector{Tuple{Int, elem_type(field)}}()
    for j in (index - 1):-1:1
        fj = polys[j]
        iszero(fj) && continue
        cj = coeff(fj, leadi)
        # If fj contains the lead of fi
        iszero(cj) && continue
        cji = div(cj, ci)
        fj = fj - cji * fi
        polys[j] = fj
        push!(multipliers, (j, -cji))
    end
    return multipliers
end

"""
    linear_relations_over_a_field(polys, preimages)

Finds all linear relations between `polys`. Encodes them in terms of the
corresponding elements of `preimages` (assuming `preimages[i] --> polys[i]`).

Returns a triple (`relations`, `polys`, `preimages`).
"""
function linear_relations_over_a_field(polys, preimages)
    @assert !isempty(polys)
    @assert length(polys) == length(preimages)
    ring = parent(first(polys))
    preimage_ring = parent(first(preimages))
    field = base_ring(ring)
    relations = Vector{elem_type(preimage_ring)}()
    # Filter out and stash zero polynomials
    permutation = collect(1:length(polys))
    zero_inds = filter(i -> iszero(polys[i]), permutation)
    for ind in zero_inds
        push!(relations, preimages[ind])
    end
    @debug "Zeroed polynomials are" preimages[zero_inds]
    permutation = setdiff(permutation, zero_inds)
    # Sort, the first monom is the smallest
    lead_monoms = map(f -> iszero(f) ? one(f) : leading_monomial(f), polys)
    sort!(permutation, by = i -> lead_monoms[i])
    polys = polys[permutation]
    preimages = preimages[permutation]
    lead_monoms = lead_monoms[permutation]
    n = length(polys)
    @inbounds for i in 1:n
        fi, multipliers = reduce_step_forward(polys, lead_monoms, i, field)
        preimage = zero(preimage_ring)
        for (idx, coef) in multipliers
            preimage += coef * preimages[idx]
        end
        preimages[i] += preimage
        polys[i] = fi
        lead_monoms[i] = iszero(fi) ? one(fi) : leading_monomial(fi)
        if iszero(fi)
            push!(relations, preimages[i])
            continue
        end
        multipliers = reduce_step_backward!(polys, lead_monoms, i, field)
        for (idx, coef) in multipliers
            preimages[idx] += coef * preimages[i]
        end
    end
    relations = map(f -> divexact(f, leading_coefficient(f)), relations)
    return relations, polys, preimages
end

function intersect_relations_over_a_field(
    relations::Vector{T},
    other_relations::Vector{T},
) where {T}
    @assert !isempty(relations)
    @assert !any(iszero, relations) && !any(iszero, other_relations)
    @assert parent(first(relations)) == parent(first(other_relations))
    ring = parent(first(relations))
    field = base_ring(ring)
    common_relations = Vector{T}()
    all_relations = vcat(relations, other_relations)
    permutation = collect(1:length(all_relations))
    lead_monoms = map(f -> leading_monomial(f), all_relations)
    sort!(permutation, by = i -> lead_monoms[i])
    all_relations = all_relations[permutation]
    lead_monoms = lead_monoms[permutation]
    n = length(all_relations)
    #
    relations_mirrored = Vector{elem_type(ring)}(undef, n)
    for i in 1:n
        absolute_index = permutation[i]
        if absolute_index <= length(relations)
            relations_mirrored[i] = all_relations[i]
        else
            relations_mirrored[i] = zero(ring)
        end
    end
    @inbounds for i in 1:n
        fi, multipliers = reduce_step_forward(all_relations, lead_monoms, i, field)
        relation_mirrored = zero(ring)
        for (j, mult) in multipliers
            relation_mirrored += mult * relations_mirrored[j]
        end
        relations_mirrored[i] += relation_mirrored
        all_relations[i] = fi
        lead_monoms[i] = iszero(fi) ? one(fi) : leading_monomial(fi)
        if iszero(fi)
            push!(common_relations, relations_mirrored[i])
            continue
        end
        multipliers = reduce_step_backward!(all_relations, lead_monoms, i, field)
        for (j, mult) in multipliers
            relations_mirrored[j] += mult * relations_mirrored[i]
        end
    end
    common_relations = map(f -> divexact(f, leading_coefficient(f)), common_relations)
    common_relations
end

"""
    linear_relations_between_normal_forms(fracs, up_to_degree)

Returns linear relations between the given `fracs` (potentially, not all).
Relations may include monomials up to the total degree `up_to_degree`.

Note: uses Monte-Carlo probabilistic algorithm. The probability of correctness
is not specified but is assumed to be close to 1.
"""
function linear_relations_between_normal_forms(
    fracs::Vector{Generic.Frac{T}},
    up_to_degree::Integer;
    seed = 42,
) where {T}
    mqs = IdealMQS(fractions_to_dennums(fracs))
    ring = parent(mqs)
    ring_param = ParamPunPam.parent_params(mqs)
    xs_param = gens(ring_param)
    nparams = nvars(ring_param)
    finite_field = Nemo.GF(2^30 + 3)
    ParamPunPam.reduce_mod_p!(mqs, finite_field)
    @info """
    Computing normal forms (probabilistic)
    Variables ($nparams in total): $xs_param
    Up to degree: $up_to_degree
    Modulo: $finite_field"""
    # We first compute relations between the normal forms of linear monomials.
    # Then, we use this knowledge to drop out some monomials of higher degrees.
    tref = TinyRowEchelonForm{Int}()
    normal_forms_ff_1, monoms_ff_1 = local_normal_forms(mqs, finite_field, 1)
    relations_ff_1 = empty(monoms_ff_1)
    for i in 1:length(normal_forms_ff_1)
        !iszero(normal_forms_ff_1[i]) && continue
        !(length(monoms_ff_1[i]) == 1) && continue
        @debug "Registering existing monomial $(monoms_ff_1[i])"
        push!(relations_ff_1, monoms_ff_1[i])
        push!(tref, exponent_vector(monoms_ff_1[i], 1))
    end
    complete_intersection_relations_ff = Vector{Nemo.gfp_mpoly}(undef, 0)
    npoints = 0
    # Compute relations at several random points until a consensus is reached
    while true
        npoints += 1
        @debug "Used specialization points: $npoints"
        @debug "Computing normal forms to to degree $up_to_degree"
        normal_forms_ff, monoms_ff =
            local_normal_forms(mqs, finite_field, up_to_degree, stop_vectors = tref)
        if isempty(normal_forms_ff)
            break
        end
        @debug "Computing relations of $(length(normal_forms_ff)) normal forms"
        relations_ff, normal_forms_ff, monoms_ff =
            linear_relations_over_a_field(normal_forms_ff, monoms_ff)
        @debug "Obtained $(length(relations_ff)) local relations"
        # do some bookkeeping to ensure that neither the accumulated
        # intersection nor the newly acquired relations are empty
        if npoints == 1
            # first point -- take all relations
            complete_intersection_relations_ff = relations_ff
            continue
        end
        if isempty(complete_intersection_relations_ff)
            # if the intersection is empty, then there are no relations
            break
        end
        if isempty(relations_ff)
            # if the newly acquired relations is empty, then there are no relations
            empty!(complete_intersection_relations_ff)
            break
        end
        @assert !isempty(complete_intersection_relations_ff)
        @assert !isempty(relations_ff)
        @debug "Intersecting relations accumulated so far" complete_intersection_relations_ff
        n_relations_ff = length(complete_intersection_relations_ff)
        complete_intersection_relations_ff = intersect_relations_over_a_field(
            complete_intersection_relations_ff,
            relations_ff,
        )
        @debug "There are $(length(complete_intersection_relations_ff)) relations in the intersection"
        m_relations_ff = length(complete_intersection_relations_ff)
        if n_relations_ff == m_relations_ff || isempty(complete_intersection_relations_ff)
            break
        end
    end
    union!(complete_intersection_relations_ff, relations_ff_1)
    @debug "Reconstructing relations to rationals"
    relations_qq = Vector{Generic.Frac{elem_type(ring_param)}}(
        undef,
        length(complete_intersection_relations_ff),
    )
    for i in 1:length(complete_intersection_relations_ff)
        relation_ff = complete_intersection_relations_ff[i]
        success, relation_qq =
            ParamPunPam.rational_reconstruct_polynomial(ring, relation_ff)
        if !success
            @warn """
            Failed to reconstruct the $i-th relation. Error will follow.
            relation: $relation_ff
            modulo: $finite_field"""
            throw(ErrorException("Rational reconstruction failed."))
        end
        relation_qq_param = evaluate(relation_qq, vcat(xs_param, one(ring)))
        relations_qq[i] = relation_qq_param // one(relation_qq_param)
    end
    relations_qq
end
