
"""
    IdealMQS

Relations between the elements of a rational functions subfield of Q(x) encoded
with an ideal in Q(x)[y]:

    <num_i(x) den_i(y) - num_i(y) den_i(x)> : Q(y)^inf

## Interface 

`IdealMQS` provides all of the functions the interface `AbstractBlackboxIdeal`
requires. See `ParamPunPam.AbstractBlackboxIdeal` for details.

## Reference

Definition 2.16
https://mediatum.ub.tum.de/doc/685465/685465.pdf
"""
mutable struct IdealMQS{T} <: AbstractBlackboxIdeal
    funcs_den_nums::Vector{Vector{T}}
    den_lcm_orig::T
    parent_ring_param::Generic.MPolyRing{T}
    # Numerators and denominators over QQ
    # NB: this lists may have different length
    nums_qq::Vector{T}
    dens_qq::Vector{T}
    sat_qq::T
    # dens_indices[i] is a pair of form (from, to).
    # Denominator as index i corresponds to numerators at indices [from..to]
    dens_indices::Vector{Tuple{Int, Int}}
    # Indices of pivot elements for each component of funcs_den_nums 
    pivots_indices::Vector{Int}
    den_lcm::T
    sat_var_index::Int
    # Numerators and denominators over GF. 
    # We cache them and maintain a map:
    # a finite field --> an image over this finite field
    cached_nums_gf::Dict{Any, Any}
    cached_dens_gf::Dict{Any, Any}
    cached_sat_gf::Dict{Any, Any}
    # Cached GBs. A mapping
    # (monomial ordering, degree) --> a GB
    cached_groebner_bases::Dict{Any, Any}

    """
        IdealMQS(funcs_den_nums::Vector{Vector})

    Given an array of polynomials, that is generators of form
        
        [[f1, f2, f3, ...], [g1, g2, g3, ...], ...]

    constructs an MQS ideal.
    """
    function IdealMQS(
        funcs_den_nums::Vector{Vector{PolyQQ}};
        sat_varname = "t",
        ordering = :degrevlex,
    ) where {PolyQQ}
        # We are given polynomials of form
        # [[f1, f2, f3, ...], [g1, g2, g3, ...], ...]
        # We prepare and store them to construct ideal specializations later 
        @assert !isempty(funcs_den_nums)
        ordering !== :degrevlex && (@warn "Ordering is not degrevlex but $ordering")
        ring = parent(first(first(funcs_den_nums)))
        @debug "Constructing the MQS ideal in $ring"
        K, n = base_ring(ring), nvars(ring)
        @debug "Finding pivot polynomials"
        # In the component f1,f2,... find the polynomial with the minimal total
        # degree. Such element will serve as a normalizing term for the
        # component
        funcs_den_nums = map(plist -> filter(!iszero, plist), funcs_den_nums)
        # funcs_den_nums = filter(plist -> length(plist) > 1, funcs_den_nums)
        @assert !isempty(funcs_den_nums) "All elements of the ideal are zero"
        pivots = map(plist -> findmin(total_degree, plist)[2], funcs_den_nums)
        pivots_indices = map(last, pivots)
        @debug "\tDegrees are $(map(first, pivots))"
        den_lcm = mapreduce(
            i -> funcs_den_nums[i][pivots_indices[i]],
            lcm,
            1:length(funcs_den_nums),
        )
        @debug "Rational functions common denominator" den_lcm
        is_constant(den_lcm) &&
            (@debug "Common denominator of the field generators is constant")
        existing_varnames = map(String, symbols(ring))
        ystrs = ["y$i" for i in 1:length(existing_varnames)]
        @assert !(sat_varname in ystrs) "The name of the saturation variable collided with a primary variable"
        sat_var_index = length(ystrs) + 1
        varnames = push!(ystrs, sat_varname)
        @debug "Saturating variable is $sat_varname, index is $sat_var_index"
        R_sat, v_sat = Nemo.PolynomialRing(K, varnames, ordering = ordering)
        # Saturation
        @assert sat_var_index == length(v_sat)
        t_sat = v_sat[sat_var_index]
        den_lcm_orig = den_lcm
        den_lcm = parent_ring_change(den_lcm, R_sat, matching = :byindex)
        den_lcm_sat = parent_ring_change(den_lcm, R_sat)
        sat_qq = den_lcm_sat * t_sat - 1
        # We construct the array of numerators nums_qq and the array of
        # denominators dens_qq. Since there are usually more unique numerators
        # than denominators, we condense the array dens_qq and produce an
        # additional array dens_indices. The element dens_indices[i] tells us
        # the indices of the numerators that correspond to a given denominator
        # dens_qq[i]
        nums_qq = empty(funcs_den_nums[1])
        dens_qq = empty(funcs_den_nums[1])
        dens_indices = Vector{Tuple{Int, Int}}()
        for i in 1:length(funcs_den_nums)
            # TODO: we can remove duplicates in numerators. Check if this helps
            plist = funcs_den_nums[i]
            den = plist[pivots_indices[i]]
            den = parent_ring_change(den, R_sat, matching = :byindex)
            push!(dens_qq, den)
            push!(dens_indices, (length(nums_qq) + 1, length(nums_qq) + length(plist) - 1))
            for j in 1:length(plist)
                j == pivots_indices[i] && continue
                num = plist[j]
                num = parent_ring_change(num, R_sat, matching = :byindex)
                push!(nums_qq, num)
            end
        end
        parent_ring_param, _ = PolynomialRing(ring, varnames, ordering = ordering)
        @debug "Constructed MQS ideal in $R_sat with $(length(nums_qq) + 1) elements"
        @assert length(pivots_indices) == length(dens_indices) == length(dens_qq)
        @assert length(pivots_indices) == length(funcs_den_nums)

        new{elem_type(R_sat)}(
            funcs_den_nums,
            den_lcm_orig,
            parent_ring_param,
            nums_qq,
            dens_qq,
            sat_qq,
            dens_indices,
            pivots_indices,
            den_lcm,
            sat_var_index,
            Dict(),
            Dict(),
            Dict(),
            Dict(),
        )
    end
end

Base.length(ideal::IdealMQS) = length(ideal.nums_qq) + 1
AbstractAlgebra.base_ring(ideal::IdealMQS) = base_ring(ideal.nums_qq[1])
AbstractAlgebra.parent(ideal::IdealMQS) = ideal.parent_ring_param
ParamPunPam.parent_params(ideal::IdealMQS) = base_ring(ideal.parent_ring_param)

@noinline function __throw_unlucky_evaluation(msg)
    throw(AssertionError("""
    Encountered a very unlucky evaluation point.
    This should not happen normally.
    (The probability of that happening is roughly 1 to 10^18).
    Please consider submitting a Github issue.

    $msg
    """))
end

# Used only for debugging!
function ideal_generators_raw(mqs::IdealMQS)
    return field_to_ideal(mqs.funcs_den_nums)
end

# Used only for debugging!
function saturate(I, Q; varname = "t", append_front = true)
    @info "Saturating the ideal, saturating variable is $varname"
    R = parent(first(I))
    K, ord = base_ring(R), ordering(R)
    existing_varnames = map(String, symbols(R))
    @assert !(varname in existing_varnames)
    varnames =
        append_front ? pushfirst!(existing_varnames, varname) :
        push!(existing_varnames, varname)
    Rt, vt = Nemo.PolynomialRing(K, varnames, ordering = ord)
    if append_front
        xs, t = vt[2:end], first(vt)
    else
        xs, t = vt[1:(end - 1)], last(vt)
    end
    It = map(f -> parent_ring_change(f, Rt), I)
    Qt = parent_ring_change(Q, Rt)
    sat = 1 - Qt * t
    push!(It, sat)
    It, t
end

# Used only for debugging!
function field_to_ideal(
    funcs_den_nums::Vector{Vector{T}};
    top_level_var = "y",
    top_level_ord = :degrevlex,
) where {T}
    @assert !isempty(funcs_den_nums)
    R = parent(first(first(funcs_den_nums)))
    @info "Producing the ideal generators in $R"
    K, n = base_ring(R), nvars(R)
    Q = reduce(lcm, map(first, funcs_den_nums))
    @debug "Rational functions common denominator" Q
    ystrs = ["$top_level_var$i" for i in 1:n]
    Ry, ys = Nemo.PolynomialRing(R, ystrs, ordering = top_level_ord)
    Qy = parent_ring_change(Q, Ry, matching = :byindex)
    I = empty(ys)
    for component in funcs_den_nums
        pivot = component[1]
        for i in 2:length(component)
            f = component[i]
            fy = parent_ring_change(f, Ry, matching = :byindex)
            qy = parent_ring_change(pivot, Ry, matching = :byindex)
            F = fy * Q - f * qy
            push!(I, F)
        end
    end
    I, t = saturate(I, Qy)
    I_rat = map(f -> map_coefficients(c -> c // one(R), f), I)
    return I_rat
end

# TODO: check that the reduction is lucky.
function ParamPunPam.reduce_mod_p!(
    mqs::IdealMQS,
    ff::Field,
) where {Field <: Union{Nemo.GaloisField, Nemo.GaloisFmpzField}}
    @debug "Reducing MQS ideal modulo $(ff)"
    # If there is a reduction modulo this field already,
    if haskey(mqs.cached_nums_gf, ff)
        @debug "Cache hit with $(ff)!"
        return nothing
    end
    nums_qq, dens_qq, sat_qq = mqs.nums_qq, mqs.dens_qq, mqs.sat_qq
    sat_gf = map_coefficients(c -> ff(c), sat_qq)
    ring_ff = parent(sat_gf)
    nums_gf = map(poly -> map_coefficients(c -> ff(c), poly, parent = ring_ff), nums_qq)
    dens_gf = map(poly -> map_coefficients(c -> ff(c), poly, parent = ring_ff), dens_qq)
    mqs.cached_nums_gf[ff] = nums_gf
    mqs.cached_dens_gf[ff] = dens_gf
    mqs.cached_sat_gf[ff] = sat_gf
    return nothing
end

function ParamPunPam.specialize_mod_p(
    mqs::IdealMQS,
    point::Vector{T};
    saturated = true,
) where {T <: Union{gfp_elem, gfp_fmpz_elem}}
    K_1 = parent(first(point))
    @debug "Evaluating MQS ideal over $K_1 at $point"
    @assert haskey(mqs.cached_nums_gf, K_1)
    nums_gf, dens_gf, sat_gf =
        mqs.cached_nums_gf[K_1], mqs.cached_dens_gf[K_1], mqs.cached_sat_gf[K_1]
    dens_indices = mqs.dens_indices
    K_2 = base_ring(nums_gf[1])
    @assert K_1 == K_2
    @assert length(point) == nvars(ParamPunPam.parent_params(mqs))
    # +1 actual variable because of the saturation!
    @assert length(point) + 1 == nvars(parent(nums_gf[1]))
    # NOTE: Assuming the saturating variable is the last one
    @assert mqs.sat_var_index == length(point) + 1
    point_sat = vcat(point, one(K_1))
    nums_gf_spec = map(num -> evaluate(num, point_sat), nums_gf)
    dens_gf_spec = map(den -> evaluate(den, point_sat), dens_gf)
    polys = Vector{typeof(sat_gf)}(undef, length(nums_gf_spec) + 1)
    @inbounds for i in 1:length(dens_gf_spec)
        den, den_spec = dens_gf[i], dens_gf_spec[i]
        iszero(den_spec) && __throw_unlucky_evaluation("Ideal: $mqs\nPoint: $point")
        span = dens_indices[i]
        for j in span[1]:span[2]
            num, num_spec = nums_gf[j], nums_gf_spec[j]
            polys[j] = num * den_spec - den * num_spec
        end
    end
    polys[end] = sat_gf
    if !saturated
        resize!(polys, length(polys) - 1)
    end
    return polys
end

function specialize(mqs::IdealMQS, point::Vector{Nemo.fmpq}; saturated = true)
    @debug "Evaluating MQS ideal over QQ at $point"
    nums_qq, dens_qq, sat_qq = mqs.nums_qq, mqs.dens_qq, mqs.sat_qq
    dens_indices = mqs.dens_indices
    K = base_ring(mqs)
    @assert length(point) == nvars(ParamPunPam.parent_params(mqs))
    # +1 actual variable because of the saturation!
    @assert length(point) + 1 == nvars(parent(nums_qq[1]))
    # NOTE: Assuming the saturating variable is the last one
    @assert mqs.sat_var_index == length(point) + 1
    point_sat = vcat(point, one(K))
    nums_qq_spec = map(num -> evaluate(num, point_sat), nums_qq)
    dens_qq_spec = map(den -> evaluate(den, point_sat), dens_qq)
    polys = Vector{typeof(sat_qq)}(undef, length(nums_qq_spec) + 1)
    @inbounds for i in 1:length(dens_qq_spec)
        den, den_spec = dens_qq[i], dens_qq_spec[i]
        iszero(den_spec) && __throw_unlucky_evaluation("Ideal: $mqs\nPoint: $point")
        span = dens_indices[i]
        for j in span[1]:span[2]
            num, num_spec = nums_qq[j], nums_qq_spec[j]
            polys[j] = num * den_spec - den * num_spec
        end
    end
    polys[end] = sat_qq
    if !saturated
        resize!(polys, length(polys) - 1)
    end
    return polys
end
