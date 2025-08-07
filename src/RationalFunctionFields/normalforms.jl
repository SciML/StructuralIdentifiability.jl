# abstract function to preform Gaussian elimination

# reducer is assumed to have the pivot coordinate equal to 1
function abstract_reduce!(prey, reducer)
    c = coordinate(prey, pivot(reducer))
    if !iszero(c)
        elementary_operation!(prey, reducer, c)
    end
end

function abstract_rref(vects::Vector)
    result = empty(vects)
    for v in vects
        for u in result
            abstract_reduce!(v, u)
        end
        if !iszero(v)
            scale!(v, 1 // coordinate(v, pivot(v)))
            for u in result
                abstract_reduce!(u, v)
            end
            push!(result, v)
        end
    end
    return result
end

# ------------------------------------------------------------------------------

# a structure with a vector of polynomials implementing the interface to fit into the above elimination algorithm

mutable struct PolyVect
    polys::Vector{MPolyRingElem}
end

function Base.iszero(p::PolyVect)
    return all(map(iszero, p.polys))
end

function pivot(p::PolyVect)
    ind = findfirst(!iszero, p.polys)
    return (ind, leading_monomial(p.polys[ind]))
end

function coordinate(p::PolyVect, index::Tuple{Int, <:MPolyRingElem})
    @assert 0 < index[1] && index[1] <= length(p.polys)
    return coeff(p.polys[index[1]], index[2])
end

function scale!(p::PolyVect, c)
    p.polys .*= c
end

function elementary_operation!(u::PolyVect, v::PolyVect, c)
    u.polys .-= (c .* v.polys)
end

# ------------------------------------------------------------------------------

function refine_relations(
    relations::Vector{T},
    gbs::Vector{Vector{T}},
) where {T <: MPolyRingElem}
    n = length(gbs)
    nf_poly_pairs = Vector{PolyVect}()
    for p in relations
        nfs = empty(relations)
        for gb in gbs
            _, nf = divrem(p, gb)
            push!(nfs, nf - constant_coefficient(nf))
        end
        push!(nfs, p)
        push!(nf_poly_pairs, PolyVect(nfs))
    end
    rref_time = @elapsed result = Vector{T}(
        map(
            pp -> pp.polys[end],
            filter(
                pp -> all(map(iszero, pp.polys[1:(end - 1)])),
                abstract_rref(nf_poly_pairs),
            ),
        ),
    )
    @debug "Rref time: $rref_time"
    return result
end

"""
    polynomial_generators(rff, up_to_degree)

Returns polynomials of degree not exceeding `up_to_degree` belonging to
the field `rff`. These polynomials generate (with products and linear combinations) all
the polynomials in the field of the degree at most `up_to_degree`.

Note: uses Monte-Carlo probabilistic algorithm. The probability of correctness
is not specified but is assumed to be close to 1.
"""
@timeit _to function polynomial_generators(
    rff::RationalFunctionField,
    up_to_degree::Integer;
    seed = 42,
)
    time_start = time_ns()
    gb_time = 0
    mqs = rff.mqs
    finite_field = Nemo.Native.GF(2^30 + 3)
    ParamPunPam.reduce_mod_p!(mqs, finite_field)

    ring_ff = parent(first(mqs.cached_nums_gf[finite_field]))
    xs_ff = contract_point(gens(ring_ff), mqs)

    all_monomials = empty(xs_ff)
    for deg in 1:up_to_degree
        for combination in Combinatorics.with_replacement_combinations(xs_ff, deg)
            push!(all_monomials, prod(combination))
        end
    end

    relations = empty(xs_ff)
    monomials_to_skip = Set([one(ring_ff)])
    for deg in 1:up_to_degree
        slice_relations =
            [m for m in setdiff(all_monomials, monomials_to_skip) if total_degree(m) <= deg]
        @debug "At degree $deg we start with $(length(slice_relations)) monomials"
        @debug "Totally $(length([m for m in all_monomials if total_degree(m) <= deg]))"
        @debug "To skip $(length([m for m in monomials_to_skip if total_degree(m) <= deg]))"

        prev_dim = 0
        curr_dim = length(slice_relations)
        batch_size = min(20, curr_dim)
        npoints = 0

        while prev_dim != curr_dim
            gbs = empty([empty(xs_ff)])
            for _ in 1:batch_size
                point = ParamPunPam.distinct_nonzero_points(finite_field, length(xs_ff))
                point_ext = extend_point(point, mqs)
                gens_spec = specialize_mod_p(mqs, point)
                gb_time += @elapsed gb_spec = Groebner.groebner(gens_spec)
                npoints += 1
                push!(gbs, gb_spec)
            end
            slice_relations = refine_relations(slice_relations, gbs)
            prev_dim, curr_dim = curr_dim, length(slice_relations)
            @debug "Current dim: $curr_dim"
        end

        append!(relations, slice_relations)
        @debug "Number of polynomials in degree $deg: $(length(slice_relations))"
        for rel in slice_relations
            new_monomials_to_skip = copy(monomials_to_skip)
            lm = leading_monomial(rel)
            for m in monomials_to_skip
                exp = 1
                while total_degree(lm) * exp + total_degree(m) <= up_to_degree
                    push!(new_monomials_to_skip, m * lm^exp)
                    exp += 1
                end
            end
            monomials_to_skip = new_monomials_to_skip
        end
        @debug "Number of points used for degree $deg: $npoints"
        _runtime_logger[:id_npoints_normalform] += npoints
    end

    @debug "Reconstructing relations to rationals"
    ring_param = ParamPunPam.parent_params(mqs)
    to_param = extend_point(gens(ring_param), mqs)
    relations_qq =
        Vector{Generic.FracFieldElem{elem_type(ring_param)}}(undef, length(relations))
    for i in 1:length(relations)
        relation_ff = relations[i]
        success, relation_qq =
            ParamPunPam.rational_reconstruct_polynomial(parent(mqs), relation_ff)
        if !success
            @warn """
            Failed to reconstruct the $i-th relation. Error will follow.
            relation: $relation_ff
            modulo: $finite_field"""
            throw(ErrorException("Rational reconstruction failed."))
        end
        relation_qq_param = evaluate(relation_qq, to_param)
        relations_qq[i] = relation_qq_param // one(relation_qq_param)
    end
    @debug "Found $(length(relations_qq)) polynomial generators using $(_runtime_logger[:id_npoints_normalform]) specialization"
    @debug "Time spent on GB computation: $gb_time"
    @info "Search for polynomial generators concluded in $((time_ns() - time_start) / 1e9)"
    _runtime_logger[:id_normalforms_time] += (time_ns() - time_start) / 1e9
    return relations_qq
end

# ------------------------------------------------------------------------------
