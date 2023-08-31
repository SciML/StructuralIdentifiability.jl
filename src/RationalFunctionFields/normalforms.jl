
"""
    linear_relations_between_normal_forms(rff, up_to_degree)

Returns the generators of the rational function field `rff` obtained as
relations over the rationals between the normal forms of the monomials up to the
degree.
"""
function linear_relations_between_normal_forms(
    rff::RationalFunctionField{T},
    up_to_degree::Integer;
    seed = 42,
) where {T}
    @assert up_to_degree > 0
    time_start = time_ns()
    # NOTE: this is not fair regarding mutation and `!`
    groebner_basis_coeffs(rff)
    gb = first(values(rff.mqs.groebner_bases))
    R = parent(gb[1])
    R_param = base_ring(base_ring(R))
    xs = gens(R)
    xs_param = gens(R_param)
    # TODO: A dirty hack!
    @assert rff.mqs.sat_var_index == length(xs)
    xs = xs[1:(end - 1)]
    @info "Computing normal forms of monomials in $(length(xs)) variables up to degree $up_to_degree"
    normal_forms = Vector{elem_type(R)}(undef, 0)
    monoms = Vector{elem_type(R_param)}(undef, 0)
    @info "GB is" gb
    @info """
    The variables rings are:
    nf. parent = $(R)
    parametric parent = $(R_param)
    gb parent = $(parent(gb[1]))"""
    @assert R == parent(gb[1])
    @assert R_param == base_ring(base_ring(parent(gb[1])))
    for deg in 1:up_to_degree
        for combination in Combinatorics.with_replacement_combinations(xs, deg)
            monom = prod(combination)
            monom_param = evaluate(monom, vcat(xs_param, one(R_param)))
            monom_mqs = monom - monom_param
            @info "Computing the normal form of" monom_mqs
            b, nf = divrem(monom_mqs, gb)
            @info "The normal form is" nf, b
            push!(monoms, numerator(monom_param))
            push!(normal_forms, nf)
        end
    end
    for a in zip(monoms, normal_forms)
        @info a
    end
    @info "Reducing the normal forms of $(length(monoms)) monomials over QQ"
    relations, normal_forms, monoms = relations_over_qq(normal_forms, monoms)
    _runtime_logger[:id_normalforms_time] = (time_ns() - time_start) / 1e9
    @info "Relations from normal forms" relations
    relations, normal_forms, monoms
end

function linear_relations_between_normal_forms_mod_p(
    rff::RationalFunctionField{T},
    up_to_degree::Integer;
    seed = 42,
) where {T}
    @assert up_to_degree > 0
    time_start = time_ns()
    groebner_basis_coeffs(rff)
    gb = first(values(rff.mqs.groebner_bases))
    basis_coeffs = map(collect ∘ coefficients, gb)
    fracs = collect(mapreduce(Set, union!, basis_coeffs))
    new_rff = RationalFunctionField(fracs)
    field = Nemo.GF(2^62 + 135)
    ParamPunPam.reduce_mod_p!(new_rff.mqs, field)
    ParamPunPam.specialize_mod_p(new_rff.mqs)
    R = parent(gb[1])
    R_param = base_ring(base_ring(R))
    xs = gens(R)
    xs_param = gens(R_param)
    # TODO: A dirty hack!
    @assert rff.mqs.sat_var_index == length(xs)
    xs = xs[1:(end - 1)]
    @info "Computing normal forms of monomials in $(length(xs)) variables up to degree $up_to_degree"
    normal_forms = Vector{elem_type(R)}(undef, 0)
    monoms = Vector{elem_type(R_param)}(undef, 0)
    @info "GB is" gb
    @info """
    The variables rings are:
    nf. parent = $(R)
    parametric parent = $(R_param)
    gb parent = $(parent(gb[1]))"""
    @assert R == parent(gb[1])
    @assert R_param == base_ring(base_ring(parent(gb[1])))
    for deg in 1:up_to_degree
        for combination in Combinatorics.with_replacement_combinations(xs, deg)
            monom = prod(combination)
            monom_param = evaluate(monom, vcat(xs_param, one(R_param)))
            monom_mqs = monom - monom_param
            @info "Computing the normal form of" monom_mqs
            b, nf = divrem(monom_mqs, gb)
            @info "The normal form is" nf, b
            push!(monoms, numerator(monom_param))
            push!(normal_forms, nf)
        end
    end
    for a in zip(monoms, normal_forms)
        @info a
    end
    @info "Reducing the normal forms of $(length(monoms)) monomials over QQ"
    relations, normal_forms, monoms = relations_over_qq(normal_forms, monoms)
    _runtime_logger[:id_normalforms_time] = (time_ns() - time_start) / 1e9
    @info "Relations from normal forms" relations
    relations, normal_forms, monoms
end

function relations_over_qq(polys, preimages)
    @assert !isempty(polys)
    fracfield = base_ring(first(polys))
    qq_relations = Vector{elem_type(fracfield)}()
    # Filter out and stash zero polynomials
    permutation = collect(1:length(polys))
    zero_inds = filter(i -> iszero(polys[i]), permutation)
    for ind in zero_inds
        push!(qq_relations, fracfield(preimages[ind]))
    end
    @debug "Zeroed monomials are" preimages[zero_inds]
    permutation = setdiff(permutation, zero_inds)
    # Sort, the first monom is the smallest
    sort!(permutation, by = i -> leading_monomial(polys[i]))
    polys = polys[permutation]
    preimages = preimages[permutation]
    lead_monoms = map(leading_monomial, polys)
    n = length(polys)
    # Polynomials live in QQ(params)[vars].
    # The first several elements are de facto elements of QQ(params).
    # NOTE: `coeff(f, i)` of a polynomial f in QQ(a)[x] is excruciatingly slow
    @inbounds for i in 1:n
        fi = polys[i]
        @debug "Reducing $i-th polynomial over QQ" fi
        qq_multipliers = map(_ -> zero(Nemo.QQ), 1:n)
        qq_multipliers[i] = one(Nemo.QQ)
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
            # If the result of division belongs to QQ.
            !is_rational_func_const(cij) && continue
            @debug "reducing $fi with $cij x $fj"
            fi = fi - cij * fj
            qq_multipliers[j] = -coeff(numerator(cij), 1)
        end
        if iszero(fi)
            @debug "Polynomial at index $i reduced to zero"
            preimage = zero(fracfield)
            for k in 1:i
                if !iszero(qq_multipliers[k])
                    preimage += qq_multipliers[k] * preimages[k]
                end
            end
            push!(qq_relations, preimage)
        end
    end
    qq_relations, polys, preimages
end

function relations_over_qq_fast(polys, preimages)
    @assert !isempty(polys)
    R = parent(first(polys))
    fracfield = base_ring(first(polys))
    paramring = base_ring(fracfield)
    qq_relations = Vector{elem_type(fracfield)}()
    # Filter out and stash zero polynomials
    permutation = collect(1:length(polys))
    zero_inds = filter(i -> iszero(polys[i]), permutation)
    for ind in zero_inds
        push!(qq_relations, fracfield(preimages[ind]))
    end
    permutation = setdiff(permutation, zero_inds)
    polys = polys[permutation]
    preimages = preimages[permutation]
    point = [Nemo.QQ(rand(1:3)) for _ in 1:nvars(paramring)]
    R_eval, _ = PolynomialRing(Nemo.QQ, symbols(R), ordering = ordering(R))
    polys_eval = Vector{elem_type(R_eval)}(undef, length(polys))
    @inbounds for i in 1:length(polys)
        f = polys[i]
        cfs = collect(coefficients(f))
        cfs_eval = Vector{Nemo.fmpq}(undef, length(cfs))
        for j in 1:length(cfs)
            num, den = unpack_fraction(cfs[j])
            num_eval = evaluate(num, point)
            den_eval = evaluate(den, point)
            # NOTE: here there is no much trouble with unlucky cancellations
            if iszero(den_eval)
                den_eval = one(den_eval)
            end
            cfs_eval[j] = num_eval // den_eval
        end
        polys_eval[i] = R_eval(cfs_eval, collect(exponent_vectors(f)))
    end
    n = length(polys_eval)
    lead_monoms = map(poly -> iszero(poly) ? one(poly) : leading_monomial(poly), polys_eval)
    # Sort, the first monom is the smallest
    permutation = collect(1:n)
    sort!(permutation, by = i -> lead_monoms[i])
    polys_eval = polys_eval[permutation]
    preimages = preimages[permutation]
    lead_monoms = lead_monoms[permutation]
    qq_multipliers = map(_ -> zero(Nemo.QQ), 1:n)
    @inbounds for i in 1:n
        fi = polys_eval[i]
        @debug "Reducing $i-th polynomial over QQ" fi
        for j in 1:(n - 1)
            qq_multipliers[j] = zero(Nemo.QQ)
        end
        qq_multipliers[i] = one(Nemo.QQ)
        for j in (i - 1):-1:1
            iszero(fi) && break
            fj = polys_eval[j]
            iszero(fj) && continue
            leadj = lead_monoms[j]
            ci = coeff(fi, leadj)
            # If fi contains the lead of fj
            iszero(ci) && continue
            cj = leading_coefficient(fj)
            cij = div(ci, cj)
            @debug "reducing $fi with $cij x $fj"
            fi = fi - cij * fj
            qq_multipliers[j] = -cij
        end
        if iszero(fi)
            @debug "Polynomial at index $i reduced to zero"
            true_relation = zero(R)
            for k in 1:i
                if !iszero(qq_multipliers[k])
                    true_relation += qq_multipliers[k] * polys[k]
                end
            end
            if !iszero(true_relation)
                continue
            end
            preimage = zero(fracfield)
            for k in 1:i
                if !iszero(qq_multipliers[k])
                    preimage += qq_multipliers[k] * preimages[k]
                end
            end
            push!(qq_relations, preimage)
        end
    end
    qq_relations, polys, preimages
end

function relations_over_qq_probababilistic(polys, preimages)
    @assert !isempty(polys)
    ring_param = parent(first(preimages))
    qq_relations = Vector{elem_type(ring_param)}()
    # Filter out and stash zero polynomials
    permutation = collect(1:length(polys))
    zero_inds = filter(i -> all(iszero, polys[i]), permutation)
    for ind in zero_inds
        push!(qq_relations, preimages[ind])
    end
    @debug "Zeroed monomials are $(preimages[zero_inds])"
    permutation = setdiff(permutation, zero_inds)
    # Sort, the first monom is the smallest
    lead_monoms =
        map(poly -> iszero(poly[1]) ? one(poly[1]) : leading_monomial(poly[1]), polys)
    sort!(permutation, by = i -> lead_monoms[i])
    polys = polys[permutation]
    preimages = preimages[permutation]
    lead_monoms = lead_monoms[permutation]
    # TODO:
    polys = map(f -> map(ff -> divexact(ff, leading_coefficient(ff)), f), polys)
    # 
    n = length(polys)
    qq_multipliers = map(_ -> zero(Nemo.QQ), 1:n)
    # Gaussian elimitation
    @inbounds for i in 1:n
        fi = polys[i][1]
        @debug "Reducing $i-th polynomial over QQ" fi
        for j in 1:(n - 1)
            qq_multipliers[j] = zero(Nemo.QQ)
        end
        qq_multipliers[i] = one(Nemo.QQ)
        for j in (i - 1):-1:1
            iszero(fi) && break
            fj = polys[j][1]
            iszero(fj) && continue
            leadj = lead_monoms[j]
            ci = coeff(fi, leadj)
            # If fi contains the lead of fj
            iszero(ci) && continue
            cj = leading_coefficient(fj)
            cij = div(ci, cj)
            @debug "reducing $fi with $cij x $fj"
            fi = fi - cij * fj
            qq_multipliers[j] = -cij
        end
        if iszero(fi)
            @debug "Polynomial at index $i reduced to zero"
            all_specializations_reduce_to_zero = true
            for j in 2:length(polys[i])
                fij = polys[i][j]
                comb = zero(fij)
                for k in 1:i
                    if !iszero(qq_multipliers[k])
                        comb += qq_multipliers[k] * polys[k][j]
                    end
                end
                if !iszero(comb)
                    all_specializations_reduce_to_zero = false
                end
            end
            if !all_specializations_reduce_to_zero
                continue
            end
            preimage = zero(ring_param)
            for k in 1:i
                if !iszero(qq_multipliers[k])
                    preimage += qq_multipliers[k] * preimages[k]
                end
            end
            push!(qq_relations, preimage)
        end
    end
    qq_relations, polys, preimages
end

function linear_relations_between_normal_forms_probabilistic(
    rff::RationalFunctionField{T},
    up_to_degree::Integer;
    seed = 42,
) where {T}
    @assert up_to_degree > 0
    time_start = time_ns()
    # Compute bases at several random points.
    # We consider a relation between normal forms correct if holds at each of
    # the considered points
    npoints = 5
    eval_bound = 5
    mqs = rff.mqs
    ring_param = ParamPunPam.parent_params(mqs)
    points = Vector{Vector{Nemo.fmpq}}(undef, npoints)
    i = 1
    while i <= npoints
        points[i] = map(_ -> Nemo.QQ(rand(1:eval_bound)), 1:nvars(ring_param))
        if points[i] in points[1:(i - 1)]
            @info "Evaluation point $(points[i]) is repeated"
            continue
        end
        i += 1
    end
    @info "Evalaution points:\n$points"
    mqs_specialized = map(point -> specialize(mqs, point), points)
    gbs_specialized = map(mqs_specialized_i -> groebner(mqs_specialized_i), mqs_specialized)
    ring_gb = parent(gbs_specialized[1][1])
    xs_gb = gens(ring_gb)
    xs_param = gens(ring_param)
    # TODO: A dirty hack!
    @assert rff.mqs.sat_var_index == length(xs_gb)
    xs_gb = xs_gb[1:(end - 1)]
    @info """
    Computing normal forms of monomials in $(length(xs_gb)) variables,
    Up to degree $up_to_degree"""
    normal_forms = Vector{Vector{elem_type(ring_gb)}}(undef, 0)
    monoms = Vector{elem_type(ring_param)}(undef, 0)
    @info """
    GBs are:
    $(join(map(string, collect(zip(points, gbs_specialized))), "\n"))
    """
    @info """
    The polynomial rings are:
    normal forms parent = $(ring_gb)
    parameters parent = $(ring_param)
    gb parent = $(parent(gbs_specialized[1][1]))"""
    # @assert R == parent(gb[1])
    # @assert R_param == base_ring(base_ring(parent(gb[1])))
    for deg in 1:up_to_degree
        for combination in Combinatorics.with_replacement_combinations(xs_gb, deg)
            monom = prod(combination)
            monom_param = evaluate(monom, vcat(xs_param, one(ring_param)))
            push!(monoms, monom_param)
            push!(normal_forms, Vector{typeof(monom)}(undef, npoints))
            for k in 1:npoints
                point = points[k]
                monom_specialized = evaluate(monom_param, point)
                monom_mqs = monom - monom_specialized
                @debug "Computing the normal form of $monom_mqs"
                _, nf = divrem(monom_mqs, gbs_specialized[k])
                @debug "The normal form is $nf"
                normal_forms[end][k] = nf
            end
        end
    end
    @info "Reducing the normal forms of $(length(monoms)) monomials over QQ"
    relations, normal_forms, monoms =
        relations_over_qq_probababilistic(normal_forms, monoms)
    _runtime_logger[:id_normalforms_time] = (time_ns() - time_start) / 1e9
    relations = map(f -> f // one(f), relations)
    @info "Relations from normal forms" relations
    relations, normal_forms, monoms
end
