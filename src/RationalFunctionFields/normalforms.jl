
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
        flag, divisor = divides(lead_vect, lead_rref)
        !flag && return false
    end
    return true
end

# ------------------------------------------------------------------------------

function normal_forms_up_to_degree_ff(
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
    @info """
    variables in finite field: $(xs_ff)
    gb parent: $(ring_ff)
    specialized gb: $(gb_ff_spec)"""
    for deg in 1:up_to_degree
        for combination in Combinatorics.with_replacement_combinations(pivot_vectors, deg)
            exp_vect = sum(combination)
            if in(stop_vectors, exp_vect)
                @info "Skipping exponent vector $exp_vect"
                continue
            end
            monom_ff = ring_ff([one(finite_field)], [exp_vect])
            monom_ff_spec = evaluate(monom_ff, point_ff_ext)
            monom_mqs_ff_spec = monom_ff - monom_ff_spec
            divisors_ff, nf_ff = divrem(monom_mqs_ff_spec, gb_ff_spec)
            @info """
            The normal form of $monom_mqs_ff_spec is:
            normalform = $nf_ff
            divisors = $divisors_ff"""
            push!(monoms_ff, monom_ff)
            push!(normal_forms_ff, nf_ff)
        end
    end
    return normal_forms_ff, monoms_ff
end

"""
    linear_relations_between_normal_forms_mod_p(fracs, up_to_degree)

Returns linear relations between the given `fracs` (potentially, not all).
Relations may include monomials up to the `up_to_degree`.

Note: uses a Monte-Carlo probabilistic algorithm. The probability of correctness
is not specified but is assumed to be close to 1.
"""
function linear_relations_between_normal_forms_mod_p(
    fracs::Vector{Generic.Frac{T}},
    up_to_degree::Integer;
    seed = 42,
) where {T}
    mqs = IdealMQS(fractions_to_dennums(fracs))
    ring = parent(mqs)
    ring_param = ParamPunPam.parent_params(mqs)
    xs_param = gens(ring_param)
    nparams = nvars(ring_param)
    ff = Nemo.GF(2^30 + 3)
    ParamPunPam.reduce_mod_p!(mqs, ff)
    @info """
    Computing normal forms (probabilistic)
    Parameters ($nparams in total): $xs_param
    Up to degree: $up_to_degree
    Modulo: $ff"""
    # We first compute relations between the normal forms of linear monomials.
    # Then, we use this knowledge to drop out some monomials of higher degrees.
    tref = TinyRowEchelonForm{Int}()
    normal_forms_ff_1, monoms_ff_1 = normal_forms_up_to_degree_ff(mqs, ff, 1)
    relations_ff_1 = empty(monoms_ff_1)
    for i in 1:length(normal_forms_ff_1)
        !iszero(normal_forms_ff_1[i]) && continue
        !(length(monoms_ff_1[i]) == 1) && continue
        @info "Registering existing monomial $(monoms_ff_1[i])"
        push!(relations_ff_1, monoms_ff_1[i])
        push!(tref, exponent_vector(monoms_ff_1[i], 1))
    end
    complete_intersection_relations_ff = Vector{Nemo.gfp_mpoly}(undef, 0)
    npoints = 0
    while true
        npoints += 1
        @info "Used specialization points: $npoints"
        normal_forms_ff, monoms_ff =
            normal_forms_up_to_degree_ff(mqs, ff, up_to_degree, stop_vectors = tref)
        @info "Computing relations of $(length(normal_forms_ff)) normal forms"
        if isempty(normal_forms_ff)
            break
        end
        relations_ff, normal_forms_ff, monoms_ff =
            relations_over_ff(normal_forms_ff, monoms_ff)
        @info "Obtained $(length(relations_ff)) local relations over FF"
        @warn "!!!"
        println("Relations:", relations_ff)
        println(filter(rel -> degree(rel, gens(parent(rel))[end - 1]) > 0, relations_ff))
        println(normal_forms_ff)
        println(monoms_ff)
        if npoints == 1
            complete_intersection_relations_ff = relations_ff
            continue
        end
        n_relations_ff = length(complete_intersection_relations_ff)
        complete_intersection_relations_ff =
            intersect(complete_intersection_relations_ff, relations_ff)
        @info "There are $(length(complete_intersection_relations_ff)) relations in the intersection"
        m_relations_ff = length(complete_intersection_relations_ff)
        if n_relations_ff == m_relations_ff
            break
        end
    end
    union!(complete_intersection_relations_ff, relations_ff_1)
    @info "Reconstructing relations to rationals"
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
            modulo: $ff"""
            throw(ErrorException("Rational reconstruction failed."))
        end
        relation_qq_param = evaluate(relation_qq, vcat(xs_param, one(ring)))
        relations_qq[i] = relation_qq_param // one(relation_qq_param)
    end
    println(relations_qq)
    relations_qq
end

function relations_over_ff(polys, preimages)
    @assert !isempty(polys)
    ring = parent(first(polys))
    preimage_ring = parent(first(preimages))
    ff = base_ring(first(polys))
    ff_relations = Vector{elem_type(preimage_ring)}()
    # Filter out and stash zero polynomials
    permutation = collect(1:length(polys))
    zero_inds = filter(i -> iszero(polys[i]), permutation)
    for ind in zero_inds
        push!(ff_relations, preimages[ind])
    end
    @info "Zeroed monomials are" preimages[zero_inds]
    permutation = setdiff(permutation, zero_inds)
    # Sort, the first monom is the smallest
    lead_monoms = map(f -> iszero(f) ? one(f) : leading_monomial(f), polys)
    sort!(permutation, by = i -> lead_monoms[i])
    polys = polys[permutation]
    preimages = preimages[permutation]
    lead_monoms = lead_monoms[permutation]
    n = length(polys)
    @inbounds for i in 1:n
        fi = polys[i]
        @info "Reducing $i-th polynomial over FF" fi
        ff_multipliers = [(i, one(ff))]
        for j in (i - 1):-1:1
            iszero(fi) && break
            fj = polys[j]
            iszero(fj) && continue
            leadj = lead_monoms[j]
            ci = coeff(fi, leadj)
            # If fi contains the lead of fj
            iszero(ci) && continue
            cj = leading_coefficient(fj)
            cij = div(ci, cj)
            @info "reducing $fi with $(-cij) x $fj"
            fi = fi - cij * fj
            push!(ff_multipliers, (j, -cij))
        end
        @info "after reduction" fi ff_multipliers
        # !iszero(fi) && continue
        # @info "Polynomial at index $i reduced to zero"
        preimage = zero(preimage_ring)
        for k in 1:length(ff_multipliers)
            idx, coef = ff_multipliers[k]
            preimage += coef * preimages[idx]
        end
        polys[i] = fi
        preimages[i] = preimage
        lead_monoms[i] = iszero(fi) ? one(fi) : leading_monomial(fi)
        if iszero(fi)
            push!(ff_relations, preimage)
            continue
        end
        @assert !iszero(fi)
        leadi = leading_monomial(fi)
        ci = leading_coefficient(fi)
        for j in (i - 1):-1:1
            fj = polys[j]
            iszero(fj) && continue
            cj = coeff(fj, leadi)
            # If fj contains the lead of fi
            iszero(cj) && continue
            cji = div(cj, ci)
            @info "Backward reducing $fj with $(-cji) x $fi"
            fj = fj - cji * fi
            preimages[j] += -cji * preimages[i]
            polys[j] = fj
        end
    end
    ff_relations = map(f -> divexact(f, leading_coefficient(f)), ff_relations)
    return ff_relations, polys, preimages
end
