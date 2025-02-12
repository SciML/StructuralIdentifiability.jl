
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
        sat_var_position = :first,
        ordering = :degrevlex,
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
        # funcs_den_nums = filter(plist -> length(plist) > 1, funcs_den_nums)
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
            plist = funcs_den_nums[i]
            den = plist[pivots_indices[i]]
            den = parent_ring_change(
                den,
                R_sat,
                matching = :byindex,
                shift = Int(sat_var_index == 1),
            )
            # push!(dens_qq, den)
            # push!(dens_indices, (length(nums_qq) + 1, length(nums_qq) + length(plist) - 1))
	    append!(dens_indices, [(_i,_i) for _i in length(nums_qq) + 1:length(nums_qq) + length(plist) - 1])
	    append!(dens_qq, [den for _i in length(nums_qq) + 1:length(nums_qq) + length(plist) - 1])
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
	 	# _num, _den = num, den
		dens_qq[length(nums_qq)+1] = _den
		push!(nums_qq, _num)
            end
        end
        parent_ring_param, _ = polynomial_ring(ring, varnames, internal_ordering = ordering)
        @debug "Constructed MQS ideal in $R_sat with $(length(nums_qq) + 1) elements"
        # @assert length(pivots_indices) == length(dens_indices) == length(dens_qq)
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

function are_generators_zero(mqs::IdealMQS)
    return all(x -> length(x) == 1, mqs.funcs_den_nums)
end

@noinline function __throw_unlucky_evaluation(msg)
    throw(AssertionError("""
    Encountered a very unlucky evaluation point.
    This should not happen normally.
    (The probability of that happening is roughly 1 to 10^18).
    Please consider submitting a Github issue.

    $msg
    """))
end

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
    @info "$(mqs.sat_var_index) $(varnames) $ring_params $(parent(mqs.sat_qq))"
    nums_qq, dens_qq, sat_qq = mqs.nums_qq, mqs.dens_qq, mqs.sat_qq
    nums_y = map(num -> parent_ring_change(num, big_ring, matching = :byindex), nums_qq)
    dens_y = map(den -> parent_ring_change(den, big_ring, matching = :byindex), dens_qq)
    sat_y = parent_ring_change(sat_qq, big_ring, matching = :byindex)
    nums_x = map(num -> parent_ring_change(num, big_ring, matching = :byname), nums_qq)
    dens_x = map(den -> parent_ring_change(den, big_ring, matching = :byname), dens_qq)
    polys = Vector{elem_type(big_ring)}(undef, length(nums_qq) + 1)
    @inbounds for i in 1:length(dens_qq)
        den_y, den_x = dens_y[i], dens_x[i]
        span = mqs.dens_indices[i]
        for j in span[1]:span[2]
            num_y, num_x = nums_y[j], nums_x[j]
            polys[j] = num_y * den_x - den_y * num_x
        end
    end
    polys[end] = sat_y
    main_var_indices = 1:(length(varnames) + 1)
    param_var_indices = (length(varnames) + 2):length(big_vars)
    return polys, main_var_indices, param_var_indices
end

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
) where {T <: Union{fpFieldElem, FpFieldElem}}
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
    point_sat = append_at_index(point, mqs.sat_var_index, one(K_1))
    nums_gf_spec = map(num -> evaluate(num, point_sat), nums_gf)
    dens_gf_spec = map(den -> evaluate(den, point_sat), dens_gf)
    polys = Vector{typeof(sat_gf)}(undef, length(nums_gf_spec) + 1)
    @inbounds for i in 1:length(dens_gf_spec)
        den, den_spec = dens_gf[i], dens_gf_spec[i]
        iszero(den_spec) && __throw_unlucky_evaluation("Point: $point")
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

function specialize(mqs::IdealMQS, point::Vector{Nemo.QQFieldElem}; saturated = true)
    @debug "Evaluating MQS ideal over QQ at $point"
    nums_qq, dens_qq, sat_qq = mqs.nums_qq, mqs.dens_qq, mqs.sat_qq
    dens_indices = mqs.dens_indices
    K = base_ring(mqs)
    @assert length(point) == nvars(ParamPunPam.parent_params(mqs))
    # +1 actual variable because of the saturation!
    @assert length(point) + 1 == nvars(parent(nums_qq[1]))
    point_sat = append_at_index(point, mqs.sat_var_index, one(K))
    nums_qq_spec = map(num -> evaluate(num, point_sat), nums_qq)
    dens_qq_spec = map(den -> evaluate(den, point_sat), dens_qq)
    polys = Vector{typeof(sat_qq)}(undef, length(nums_qq_spec) + 1)
    @inbounds for i in 1:length(dens_qq_spec)
        den, den_spec = dens_qq[i], dens_qq_spec[i]
        iszero(den_spec) && __throw_unlucky_evaluation("Point: $point")
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
