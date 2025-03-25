
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
    # Indices of pivot elements for each component of funcs_den_nums 
    pivots_indices::Vector{Int}
    parent_ring_param::Generic.MPolyRing{T}
    # Numerators and denominators over QQ
    nums_qq::Vector{T}
    dens_qq::Vector{T}
    const_polys::Vector{T} # the first one is the sat, but could be others
    den_lcm::T
    sat_var_index::Int
    # Numerators and denominators over GF. 
    # We cache them and maintain a map:
    # a finite field --> an image over this finite field
    cached_nums_gf::Dict{Any, Any}
    cached_dens_gf::Dict{Any, Any}
    cached_const_polys_gf::Dict{Any, Any}
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
        sat_var_position = :first,
        ordering = :degrevlex,
        extra_const_polys::Vector{PolyQQ}=Vector{PolyQQ}(),
    ) where {PolyQQ}
        # We are given polynomials of form
        # [[f1, f2, f3, ...], [g1, g2, g3, ...], ...]
        # We prepare and store them to construct ideal specializations later 
        @assert !isempty(funcs_den_nums)
        @assert sat_var_position in (:first, :last)
        ordering !== :degrevlex && (@warn "Ordering is not degrevlex but $ordering")
        ring = parent(first(first(funcs_den_nums)))
        @debug "Constructing the MQS ideal in $ring"
        K, n = base_ring(ring), nvars(ring)
        @debug "Finding pivot polynomials"
        # In the component f1,f2,... find the polynomial with the minimal total
        # degree and length. Such element will serve as a normalizing term for
        # the component
        funcs_den_nums = map(plist -> filter(!iszero, plist), funcs_den_nums)
        @assert !isempty(funcs_den_nums) "All elements of the ideal are zero"
        pivots =
            map(plist -> findmin(p -> (total_degree(p), length(p)), plist), funcs_den_nums)
        pivots_indices = map(last, pivots)
        for plist in funcs_den_nums
            @debug "\tDegrees in this list are $(map(total_degree, plist))"
        end
        @debug "\tDegrees and lengths are $(map(first, pivots))"
        den_lcm = mapreduce(
            i -> funcs_den_nums[i][pivots_indices[i]],
            lcm,
            1:length(funcs_den_nums),
        )
        @debug "Rational functions common denominator is of degree $(total_degree(den_lcm)) and of length $(length(den_lcm))"
        is_constant(den_lcm) &&
            (@debug "Common denominator of the field generators is constant")
        existing_varnames = map(String, symbols(ring))
        ystrs = ["y$i" for i in 1:length(existing_varnames)]
        @assert !(sat_varname in ystrs) "The name of the saturation variable collided with a primary variable"
        sat_var_index = if sat_var_position === :first
            1
        else
            @assert sat_var_position === :last
            length(ystrs) + 1
        end
        varnames = append_at_index(ystrs, sat_var_index, sat_varname)
        @debug "Saturating variable is $sat_varname, index is $sat_var_index"
        R_sat, v_sat = Nemo.polynomial_ring(K, varnames, internal_ordering = ordering)
        # Saturation
        t_sat = v_sat[sat_var_index]
        den_lcm_orig = den_lcm
        den_lcm = parent_ring_change(
            den_lcm,
            R_sat,
            matching = :byindex,
            shift = Int(sat_var_index == 1),
        )
        den_lcm_sat = parent_ring_change(den_lcm, R_sat)
        sat_qq = den_lcm_sat * t_sat - 1
        const_polys = [sat_qq]

        # We construct the array of numerators nums_qq and the array of
        # denominators dens_qq of the same length
        nums_qq = empty(funcs_den_nums[1])
        dens_qq = empty(funcs_den_nums[1])
        for i in 1:length(funcs_den_nums)
            plist = funcs_den_nums[i]
            den = plist[pivots_indices[i]]
            den = parent_ring_change(
                den,
                R_sat,
                matching = :byindex,
                shift = Int(sat_var_index == 1),
            )
            for j in 1:length(plist)
                j == pivots_indices[i] && continue
                num = plist[j]
                num = parent_ring_change(
                    num,
                    R_sat,
                    matching = :byindex,
                    shift = Int(sat_var_index == 1),
                )
                G = gcd(num, den)
                _num, _den = num / G, den / G
                push!(nums_qq, _num)
                push!(dens_qq, _den)
            end
        end
        parent_ring_param, _ = polynomial_ring(ring, varnames, internal_ordering = ordering)
        @debug "Constructed MQS ideal in $R_sat with $(length(nums_qq) + 1) elements"
        @assert length(pivots_indices) == length(funcs_den_nums)
        @assert length(nums_qq) == length(dens_qq)

        new{elem_type(R_sat)}(
            funcs_den_nums,
            den_lcm_orig,
            pivots_indices,
            parent_ring_param,
            nums_qq,
            dens_qq,
            const_polys,
            den_lcm,
            sat_var_index,
            Dict(),
            Dict(),
            Dict(),
            Dict(),
        )
    end
end

# ------------------------------------------------------------------------------

Base.length(ideal::IdealMQS) = length(ideal.nums_qq) + 1
AbstractAlgebra.base_ring(ideal::IdealMQS) = base_ring(ideal.nums_qq[1])
AbstractAlgebra.parent(ideal::IdealMQS) = ideal.parent_ring_param
ParamPunPam.parent_params(ideal::IdealMQS) = base_ring(ideal.parent_ring_param)

# ------------------------------------------------------------------------------

function are_generators_zero(mqs::IdealMQS)
    return all(x -> length(x) == 1, mqs.funcs_den_nums)
end

# ------------------------------------------------------------------------------

@noinline function __throw_unlucky_evaluation(msg)
    throw(AssertionError("""
    Encountered a very unlucky evaluation point.
    This should not happen normally.
    (The probability of that happening is roughly 1 to 10^18).
    Please consider submitting a Github issue.

    $msg
    """))
end

# ------------------------------------------------------------------------------

function fractionfree_generators_raw(mqs::IdealMQS)
    # TODO: this assumes mqs.sat_var_index is last, and thus is broken
    ring_params = ParamPunPam.parent_params(mqs)
    K = base_ring(ring_params)
    varnames = map(string, Nemo.symbols(ring_params))
    # The hope is that new variables' names would not intersect with the old ones
    old_varnames = map(i -> "y$i", 1:length(varnames))
    new_varnames = map(i -> "__var_$i", 1:(length(varnames) + 1))
    if !isempty(intersect(old_varnames, new_varnames))
        @warn "Intersection in two sets of variables! $varnames $new_varnames"
    end
    # NOTE: new variables go first!
    big_ring, big_vars =
        polynomial_ring(K, vcat(new_varnames, old_varnames), internal_ordering = :lex)
    @info "$(mqs.sat_var_index) $(varnames) $ring_params $(parent(first(mqs.const_polys)))"
    nums_qq, dens_qq, const_polys = mqs.nums_qq, mqs.dens_qq, mqs.const_polys
    nums_y = map(num -> parent_ring_change(num, big_ring, matching = :byindex), nums_qq)
    dens_y = map(den -> parent_ring_change(den, big_ring, matching = :byindex), dens_qq)
    const_polys_y = map(p -> parent_ring_change(p, big_ring, matching = :byindex), const_polys)
    nums_x = map(num -> parent_ring_change(num, big_ring, matching = :byname), nums_qq)
    dens_x = map(den -> parent_ring_change(den, big_ring, matching = :byname), dens_qq)
    polys = Vector{elem_type(big_ring)}(undef, length(nums_qq) + length(const_polys))
    @inbounds for i in 1:length(dens_qq)
        den_y, den_x = dens_y[i], dens_x[i]
        num_y, num_x = nums_y[i], nums_x[i]
        polys[i] = num_y * den_x - den_y * num_x
    end
    @inbounds for i in 1:length(const_polys)
        polys[length(nums_qq) + i] = const_polys_y[i]
    end
    main_var_indices = 1:(length(varnames) + 1)
    param_var_indices = (length(varnames) + 2):length(big_vars)
    return polys, main_var_indices, param_var_indices
end

# ------------------------------------------------------------------------------

# TODO: check that the reduction is lucky.
function ParamPunPam.reduce_mod_p!(
    mqs::IdealMQS,
    ff::Field,
) where {Field <: Union{Nemo.fpField, Nemo.FpField}}
    @debug "Reducing MQS ideal modulo $(ff)"
    # If there is a reduction modulo this field already,
    if haskey(mqs.cached_nums_gf, ff)
        @debug "Cache hit with $(ff)!"
        return nothing
    end
    nums_qq, dens_qq, const_polys = mqs.nums_qq, mqs.dens_qq, mqs.const_polys
    ring_qq = parent(first(const_polys))
    ring_ff, _ = Nemo.polynomial_ring(ff, map(var_to_str, gens(ring_qq)), internal_ordering=Nemo.internal_ordering(ring_qq))
    nums_gf = map(poly -> map_coefficients(c -> ff(c), poly, parent = ring_ff), nums_qq)
    dens_gf = map(poly -> map_coefficients(c -> ff(c), poly, parent = ring_ff), dens_qq)
    const_polys_gf = map(poly -> map_coefficients(c -> ff(c), poly, parent = ring_ff), const_polys)
    mqs.cached_nums_gf[ff] = nums_gf
    mqs.cached_dens_gf[ff] = dens_gf
    mqs.cached_const_polys_gf[ff] = const_polys_gf
    return nothing
end

# ------------------------------------------------------------------------------

"""
    Returns a vector of MQS-polynomial evaluated at point
"""
function fractions_to_mqs_specialized(
    nums::Vector{T},
    dens::Vector{T},
    point::Vector{P},
) where {T, P}
    @assert length(nums) == length(dens)
    @assert length(gens(parent(first(nums)))) == length(point)
    polys = Vector{typeof(first(nums))}(undef, length(nums))
    nums_spec = map(poly -> evaluate(poly, point), nums)
    dens_spec = map(poly -> evaluate(poly, point), dens)
    @inbounds for i in 1:length(dens_spec)
        den, den_spec = dens[i], dens_spec[i]
        iszero(den_spec) && __throw_unlucky_evaluation("Point: $point")
        num, num_spec = nums[i], nums_spec[i]
        polys[i] = num * den_spec - den * num_spec
    end
    return polys
end

# ------------------------------------------------------------------------------

function ParamPunPam.specialize_mod_p(
    mqs::IdealMQS,
    point::Vector{T},
) where {T <: Union{fpFieldElem, FpFieldElem}}
    K_1 = parent(first(point))
    @debug "Evaluating MQS ideal over $K_1 at $point"
    @assert haskey(mqs.cached_nums_gf, K_1)
    nums_gf, dens_gf, const_polys_gf =
        mqs.cached_nums_gf[K_1], mqs.cached_dens_gf[K_1], mqs.cached_const_polys_gf[K_1]
    K_2 = base_ring(nums_gf[1])
    @assert K_1 == K_2
    @assert length(point) == nvars(ParamPunPam.parent_params(mqs))
    point_sat = append_at_index(point, mqs.sat_var_index, one(K_1))
    result = fractions_to_mqs_specialized(nums_gf, dens_gf, point_sat)
    append!(result, const_polys_gf)
    return result
end

# ------------------------------------------------------------------------------

function specialize(mqs::IdealMQS, point::Vector{Nemo.QQFieldElem})
    @debug "Evaluating MQS ideal over QQ at $point"
    nums_qq, dens_qq = mqs.nums_qq, mqs.dens_qq
    K = base_ring(mqs)
    @assert length(point) == nvars(ParamPunPam.parent_params(mqs))
    point_sat = append_at_index(point, mqs.sat_var_index, one(K))
    result = fractions_to_mqs_specialized(nums_qq, dens_qq, point_sat)
    append!(result, mqs.const_polys)
    return result
end

# ------------------------------------------------------------------------------

function specialize_fracs_to_mqs(mqs::IdealMQS, fracs, point)
    ff = parent(first(point))
    point_sat = append_at_index(point, mqs.sat_var_index, one(ff))
    new_ring = parent(mqs.nums_qq[1])
    if characteristic(ff) > 0
        @assert haskey(mqs.cached_nums_gf, ff)
        new_ring = parent(first(mqs.cached_nums_gf[ff]))
    end
    num_den_pairs = map(
        pair -> map(
            p -> parent_ring_change(
                p,
                new_ring,
                matching = :byindex,
                shift = Int(mqs.sat_var_index == 1),
            ),
            pair,
        ),
        map(unpack_fraction, fracs),
    )
    nums = map(first, num_den_pairs)
    dens = map(last, num_den_pairs)
    return fractions_to_mqs_specialized(nums, dens, point_sat)
end
