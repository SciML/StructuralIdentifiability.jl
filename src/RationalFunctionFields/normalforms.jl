
# Maintains a row echelon form of a set of vectors over the integrals.
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
    local_normal_forms(mqs, field, up_to_degree, point, [stop_vectors])

Computes the normal forms of MQS-monomials modulo the MQS-ideal `mqs`
specialized at `point`. Considers monomials up to total degree `up_to_degree`
over the `field`.

Ignores any monomials whose exponent vectors are present in `stop_vectors`.

Returns a triple (`gb`, `normal_forms`, `monomials`).
"""
function local_normal_forms(
    mqs::IdealMQS,
    finite_field,
    up_to_degree::Integer,
    point;
    stop_vectors = TinyRowEchelonForm{Int}(),
)
    @assert !isempty(point)
    @assert parent(first(point)) == finite_field
    point_ff_ext = append_at_index(point, mqs.sat_var_index, one(finite_field))
    gens_ff_spec = specialize_mod_p(mqs, point)
    gb_ff_spec = Groebner.groebner(gens_ff_spec, loglevel = _groebner_loglevel[])
    ring_ff = parent(gb_ff_spec[1])
    xs_ff = gens(ring_ff)
    normal_forms_ff = Vector{elem_type(ring_ff)}(undef, 0)
    monoms_ff = Vector{elem_type(ring_ff)}(undef, 0)
    xs_ff = cut_at_index(xs_ff, mqs.sat_var_index)
    pivot_vectors = map(f -> exponent_vector(f, 1), xs_ff)
    @debug """
    variables in the finite field: $(xs_ff)
    gb parent: $(ring_ff)
    specialized gb: $(gb_ff_spec)
    Evaluation point: $point_ff_ext"""
    # Compute the normal forms of all monomials of degrees up to `up_to_degree`
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
    return gb_ff_spec, normal_forms_ff, monoms_ff
end

# Linearly reduce polys[index] w.r.t polys[1..index-1].
# Return the result and a vector of units, multipliers of reducers.
function reduce_step_forward(
    polys::Vector{T},
    lead_monoms::Vector{T},
    index::Int,
    field,
) where {T}
    @assert length(polys) == length(lead_monoms)
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

# Linearly reduce polys[1..index-1] w.r.t polys[index], inplace.
# Return a vector of units, multipliers of reducer.
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
    @debug "Zeroed polynomials are $(preimages[zero_inds])"
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

"""
    intersect_relations_over_a_field(original_relations, other_relations)

Returns the intersection of `original_relations` and `other_relations` as linear
subspaces.
"""
function intersect_relations_over_a_field(
    original_relations::Vector{T},
    other_relations::Vector{T},
) where {T}
    if isempty(original_relations) || isempty(other_relations)
        return empty(original_relations)
    end
    @assert !any(iszero, original_relations) && !any(iszero, other_relations)
    @assert parent(first(original_relations)) == parent(first(other_relations))
    ring = parent(first(original_relations))
    field = base_ring(ring)
    common_relations = Vector{T}()
    # We compute the row echelon form of all available relations.
    #
    # If matrix row reduces to zero, or, equivaletly, the combination
    # a*original_relations + b*other_relations vanishes for some a,b, then we
    # take a*original_relations as a vector in the intersection.
    all_relations = vcat(original_relations, other_relations)
    permutation = collect(1:length(all_relations))
    lead_monoms = map(f -> leading_monomial(f), all_relations)
    sort!(permutation, by = i -> lead_monoms[i])
    all_relations = all_relations[permutation]
    lead_monoms = lead_monoms[permutation]
    n = length(all_relations)
    # We basically perform the same elementary transformations on two matrices
    # simultaneously. 
    #
    # The first matrix includes ALL relations, and it governs the process of
    # reduction. The second matrix includes only the original relations, and it
    # lazily mirrors all elementary transformations.
    mirrored_relations = Vector{Dict{Int, elem_type(field)}}(undef, n)
    for i in 1:n
        # mirrored_relations[i] is a sparse vector that keeps track of
        # polynomials that contributed to reducing the polynomial at index i
        mirrored_relations[i] = Dict{Int, elem_type(field)}()
        absolute_index = permutation[i]
        if absolute_index <= length(original_relations)
            mirrored_relations[i][i] = one(field)
        end
    end
    all_relations_copy = copy(all_relations)
    @inbounds for i in 1:n
        # Reduce the i-th row with the 1..i-1 rows
        remainder, multipliers = reduce_step_forward(all_relations, lead_monoms, i, field)
        for (j, mult) in multipliers
            for (k, mult2) in mirrored_relations[j]
                new_mult = mult * mult2
                if haskey(mirrored_relations[i], k)
                    mirrored_relations[i][k] += new_mult
                else
                    mirrored_relations[i][k] = new_mult
                end
            end
        end
        all_relations[i] = remainder
        lead_monoms[i] = iszero(remainder) ? one(ring) : leading_monomial(remainder)
        # Relation found!
        if iszero(remainder)
            relation = zero(ring)
            for (index, mult) in mirrored_relations[i]
                relation += mult * all_relations_copy[index]
            end
            push!(common_relations, relation)
            continue
        end
        # Reduce the 1..i-1 rows with the i-th row
        @assert !iszero(all_relations[i])
        multipliers = reduce_step_backward!(all_relations, lead_monoms, i, field)
        for (j, mult) in multipliers
            for (k, mult2) in mirrored_relations[i]
                new_mult = mult * mult2
                if haskey(mirrored_relations[j], k)
                    mirrored_relations[j][k] += new_mult
                else
                    mirrored_relations[j][k] = new_mult
                end
            end
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
@timeit _to function linear_relations_between_normal_forms(
    fracs::Vector{Generic.Frac{T}},
    up_to_degree::Integer;
    seed = 42,
) where {T}
    time_start = time_ns()
    mqs = IdealMQS(fractions_to_dennums(fracs))
    ring = parent(mqs)
    ring_param = ParamPunPam.parent_params(mqs)
    xs_param = gens(ring_param)
    nparams = nvars(ring_param)
    finite_field = Nemo.GF(2^30 + 3)
    ParamPunPam.reduce_mod_p!(mqs, finite_field)
    @info "Computing normal forms of degree $up_to_degree in $nparams variables"
    @debug """Variables ($nparams in total): $xs_param
    Modulo: $finite_field"""
    # We first compute relations between the normal forms of linear monomials.
    # Then, we use this knowledge to drop out some monomials of higher degrees.
    tref = TinyRowEchelonForm{Int}()
    point = ParamPunPam.distinct_nonzero_points(finite_field, nvars(ring_param))
    _, normal_forms_ff_1, monoms_ff_1 = local_normal_forms(mqs, finite_field, 1, point)
    relations_ff_1 = empty(monoms_ff_1)
    for i in 1:length(normal_forms_ff_1)
        !iszero(normal_forms_ff_1[i]) && continue
        !(length(monoms_ff_1[i]) == 1) && continue
        @debug "Registering existing monomial $(monoms_ff_1[i])"
        push!(relations_ff_1, monoms_ff_1[i])
        push!(tref, exponent_vector(monoms_ff_1[i], 1))
    end
    complete_intersection_relations_ff = Vector{Nemo.gfp_mpoly}(undef, 0)
    iters = 0
    # Compute relations at several random points until a consensus is reached
    while true
        iters += 1
        point = ParamPunPam.distinct_nonzero_points(finite_field, nvars(ring_param))
        @debug "Used specialization points: $iters"
        @debug "Computing normal forms to to degree $up_to_degree"
        gb_ff, normal_forms_ff, monoms_ff =
            local_normal_forms(mqs, finite_field, up_to_degree, point, stop_vectors = tref)
        if isempty(normal_forms_ff)
            break
        end
        @debug "Computing relations of $(length(normal_forms_ff)) normal forms"
        relations_ff, normal_forms_ff, monoms_ff =
            linear_relations_over_a_field(normal_forms_ff, monoms_ff)
        @debug "Obtained $(length(relations_ff)) local relations"
        if iters == 1
            # first point in the sequence, take all relations
            complete_intersection_relations_ff = relations_ff
            continue
        end
        if isempty(relations_ff)
            # if the newly acquired relations is empty, then the intersection is
            # empty
            empty!(complete_intersection_relations_ff)
            break
        end
        if isempty(complete_intersection_relations_ff)
            # if the intersection is empty, then there are no relations
            break
        end
        @assert !isempty(complete_intersection_relations_ff)
        @assert !isempty(relations_ff)
        n_relations_ff = length(complete_intersection_relations_ff)
        # Filter out some relations from the complete intersection
        zeroed_relations_inds = Vector{Int}()
        point_ext = append_at_index(point, mqs.sat_var_index, one(finite_field))
        for i in 1:length(complete_intersection_relations_ff)
            relation = complete_intersection_relations_ff[i]
            relation_mqs = relation - evaluate(relation, point_ext)
            _, nf = divrem(relation_mqs, gb_ff)
            if iszero(nf)
                push!(zeroed_relations_inds, i)
            end
        end
        @debug """
   Relations in the previous intersection: $(length(complete_intersection_relations_ff))
   Vanished at the current point: $(length(zeroed_relations_inds))"""
        non_zeroed_relations_inds =
            setdiff(collect(1:n_relations_ff), zeroed_relations_inds)
        zeroed_relations_from_complete_intersection =
            complete_intersection_relations_ff[zeroed_relations_inds]
        non_zeroed_relations_from_complete_intersection =
            complete_intersection_relations_ff[non_zeroed_relations_inds]
        # Fairly intersect as vector subspaces 
        non_zeroed_relations_from_complete_intersection = intersect_relations_over_a_field(
            non_zeroed_relations_from_complete_intersection,
            relations_ff,
        )
        complete_intersection_relations_ff = vcat(
            non_zeroed_relations_from_complete_intersection,
            zeroed_relations_from_complete_intersection,
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
        relation_qq_param = evaluate(
            relation_qq,
            append_at_index(xs_param, mqs.sat_var_index, one(ring_param)),
        )
        relations_qq[i] = relation_qq_param // one(relation_qq_param)
    end
    @info "Used $iters specializations in $((time_ns() - time_start) / 1e9) seconds, found $(length(complete_intersection_relations_ff)) relations"
    relations_qq
end
