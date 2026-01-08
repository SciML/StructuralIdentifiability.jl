# ------------------------------------------------------------------------------

function append_at_index(vec::Vector{T}, idx::Integer, val::T) where {T}
    # NOTE: could also just use insert!
    @assert (idx == 1) || (idx == length(vec) + 1)
    return if idx == 1
        vcat(val, vec)
    else
        vcat(vec, val)
    end
end

function cut_at_index(vec::Vector{T}, idx::Integer) where {T}
    @assert (idx == 1) || (idx == length(vec))
    return if idx == 1
        vec[2:end]
    else
        vec[1:(end - 1)]
    end
end

# ------------------------------------------------------------------------------

function simplify_frac(numer::P, denom::P) where {P <: MPolyRingElem}
    gcd_sub = gcd(numer, denom)
    sub_numer = divexact(numer, gcd_sub)
    sub_denom = divexact(denom, gcd_sub)
    return sub_numer, sub_denom
end

# ------------------------------------------------------------------------------

"""
    make_substitution(f, var_sub, val_numer, val_denom)

Substitute a variable in a polynomial with an expression

Input:
- `f` - the polynomial
- `var_sub` - the variable to be substituted
- `var_numer` - numerator of the substitution expression
- `var_denom` - denominator of the substitution expression

Output:
- `polynomial` - result of substitution
"""
function make_substitution(
        f::P,
        var_sub::P,
        val_numer::P,
        val_denom::P,
    ) where {P <: MPolyRingElem}
    d = Nemo.degree(f, var_sub)

    result = 0
    @debug "Substitution in a polynomial of degree $d"
    for i in 0:d
        @debug "\t Degree $i"
        result += coeff(f, [var_sub], [i]) * (val_numer^i) * (val_denom^(d - i))
        @debug "\t Intermediate result of size $(length(result))"
    end
    return result
end

# ------------------------------------------------------------------------------

function homogenize(fs)
    ring = parent(fs[1])
    newring, homogeneous_vars = polynomial_ring(
        base_ring(ring),
        vcat("X0", map(string, gens(ring))),
        internal_ordering = internal_ordering(ring),
    )
    Fs = empty(fs)
    for f in fs
        D = total_degree(f)
        new_f = zero(newring)
        for term in terms(f)
            cf = coeff(term, 1)
            ev = monomial(term, 1)
            d = total_degree(ev)
            new_f +=
                cf * evaluate(ev, homogeneous_vars[2:end]) * homogeneous_vars[1]^(D - d)
        end
        push!(Fs, new_f)
    end
    return Fs
end

function dehomogenize(Fs)
    ring = parent(Fs[1])
    newring, dehom_vars = polynomial_ring(
        base_ring(ring),
        map(string, gens(ring)[2:end]),
        internal_ordering = internal_ordering(ring),
    )
    fs = empty(Fs)
    for F in Fs
        f = evaluate(F, vcat(one(newring), dehom_vars))
        push!(fs, f)
    end
    return fs
end

# ------------------------------------------------------------------------------

"""
    uncertain_factorization(f)

Input:
- `f` - polynomial with rational coefficients

Output:
- list of pairs `(div, certainty)` where
    - `div`'s are divisors of `f` such that `f` is their product with certain powers
    - if `certainty` is true, `div` is ``Q``-irreducible
"""
function uncertain_factorization(f::MPolyRingElem{QQFieldElem})
    vars_f = vars(f)
    if isempty(vars_f)
        return Array{Tuple{typeof(f), Bool}, 1}()
    end
    main_var = vars_f[end]
    d = Nemo.degree(f, main_var)
    mainvar_coeffs = [coeff(f, [main_var], [i]) for i in 0:d]
    gcd_coef = mainvar_coeffs[end]
    for i in d:-1:1
        gcd_coef = gcd(gcd_coef, mainvar_coeffs[i])
    end
    f = divexact(f, gcd_coef)

    is_irr = false
    while true
        plugin = rand(-3:3, length(gens(parent(f))))
        if evaluate(mainvar_coeffs[end], plugin) != 0
            uni_ring, T = Nemo.polynomial_ring(base_ring(f), "T")
            f_uni = sum([evaluate(mainvar_coeffs[i + 1], plugin) * T^i for i in 0:d])
            if !isone(gcd(f_uni, derivative(f_uni)))
                factor_out = gcd(f, derivative(f, main_var))
                if degree(factor_out, main_var) != degree(gcd(f_uni, derivative(f_uni)))
                    continue
                end
                @debug "Nonsquarefree poly, dividing by $factor_out"
                f = divexact(f, factor_out)
                f_uni = divexact(f_uni, gcd(f_uni, derivative(f_uni)))
            end
            is_irr = Nemo.is_irreducible(f_uni)
            break
        end
    end

    coeff_factors = uncertain_factorization(gcd_coef)
    push!(coeff_factors, (f, is_irr))
    return coeff_factors
end

# ------------------------------------------------------------------------------

function fast_factor(poly::MPolyRingElem{QQFieldElem})
    prelim_factors = uncertain_factorization(poly)
    cert_factors = map(pair -> pair[1], filter(f -> f[2], prelim_factors))
    uncert_factors = map(pair -> pair[1], filter(f -> !f[2], prelim_factors))
    for p in uncert_factors
        for f in Nemo.factor(p)
            push!(cert_factors, f[1])
        end
    end
    return cert_factors
end

# ------------------------------------------------------------------------------

function dict_to_poly(dict_monom::Dict{Array{Int, 1}, <:RingElem}, poly_ring::MPolyRing)
    builder = MPolyBuildCtx(poly_ring)
    for (monom, coef) in pairs(dict_monom)
        push_term!(builder, base_ring(poly_ring)(coef), monom)
    end
    return finish(builder)
end

# ------------------------------------------------------------------------------

"""
    extract_coefficients(poly, variables)

Input:
- `poly` - multivariate polynomial
- `variables` - a list of variables from the generators of the ring of `poly`
Output:
- dictionary with keys being tuples of length `length(variables)` and values being polynomials in the variables other than those which are the coefficients at the corresponding monomials (in a smaller polynomial ring)
"""
function extract_coefficients(poly::P, variables::Array{P, 1}) where {P <: MPolyRingElem}
    xs = gens(parent(poly))
    @assert all(in(xs), variables)
    # Use a type-stable version by converting to Int explicitly
    cut_indices = Vector{Int}(undef, length(variables))
    for (j, v) in enumerate(variables)
        idx = findfirst(x -> x == v, xs)
        cut_indices[j] = idx::Int  # Assert non-nothing for type stability
    end
    coeff_indices = setdiff(collect(1:length(xs)), cut_indices)
    coeff_vars = xs[coeff_indices]

    K = base_ring(parent(poly))
    new_ring, _ = Nemo.polynomial_ring(K, map(vv -> var_to_str(vv, xs = xs), coeff_vars))
    FieldType = elem_type(K)

    result = Dict{Vector{Int}, Tuple{Vector{Vector{Int}}, Vector{FieldType}}}()

    n_cut = length(cut_indices)
    n_coeff = length(coeff_indices)
    @inbounds for i in 1:length(poly)
        coef = coeff(poly, i)
        evec = exponent_vector(poly, i)
        var_slice = Vector{Int}(undef, n_cut)
        for j in 1:n_cut
            var_slice[j] = evec[cut_indices[j]]
        end
        if !haskey(result, var_slice)
            monom_vect, coef_vect = Vector{Vector{Int}}(), Vector{FieldType}()
            sizehint!(monom_vect, 8)
            sizehint!(coef_vect, 8)
            result[var_slice] = (monom_vect, coef_vect)
        end
        monom_vect, coef_vect = result[var_slice]
        new_monom = Vector{Int}(undef, n_coeff)
        for j in 1:n_coeff
            new_monom[j] = evec[coeff_indices[j]]
        end
        push!(monom_vect, new_monom)
        push!(coef_vect, coef)
    end

    return Dict(k => new_ring(v[2], v[1]) for (k, v) in result)
end

# ------------------------------------------------------------------------------

"""
    switch_ring(v, ring)

For a variable `v`, returns a variable in `ring` with the same name
"""
function switch_ring(v::MPolyRingElem, ring::MPolyRing)
    ind = findfirst(vv -> vv == v, gens(parent(v)))

    return str_to_var(string(symbols(parent(v))[ind]), ring)
end

# ------------------------------------------------------------------------------

"""
    decompose_derivative(varname, prefixes)

Determines if it is possible to represent the `varname` as `a_number` where `a` is an element of `prefixes`
If yes, returns a pair `(a, number)`, otherwise nothing
"""
function decompose_derivative(varname::String, prefixes::Array{String})
    for pr in prefixes
        if startswith(varname, pr) && length(varname) > length(pr) + 1
            if varname[nextind(varname, ncodeunits(pr))] == '_' &&
                    all(isdigit, varname[(nextind(varname, ncodeunits(pr)) + 1):end])
                return (pr, parse(Int, varname[(nextind(varname, ncodeunits(pr)) + 1):end]))
            end
        end
    end
    return
end

# -----------------------------------------------------------------------------

"""
    difforder(diffpoly, prefix)

Finds the order of a differential polynomial `diffpoly` in a variable `prefix`,
returns -1 is the variable does not appear
"""

function difforder(diffpoly::MPolyRingElem, prefix::String)
    orders = [-1]
    for v in vars(diffpoly)
        d = decompose_derivative(var_to_str(v), [prefix])
        if d != nothing
            push!(orders, d[2])
        end
    end
    return max(orders...)
end

# -----------------------------------------------------------------------------

"""
    replace_with_ic(ode::ODE, funcs)

Takes an ode and a list of functions in the states and parameters and makes a change of variable
names `x(t) -> x(0)`. Function is used to prepare the output for the case of known initial conditions
"""
function replace_with_ic(ode, funcs)
    varnames = [(var_to_str(p), var_to_str(p)) for p in ode.parameters]
    for x in ode.x_vars
        s = var_to_str(x)
        if endswith(s, "(t)")
            push!(varnames, (s, s[1:(end - 3)] * "(0)"))
        else
            push!(varnames, (s, s * "(0)"))
        end
    end
    R0, vars0 = polynomial_ring(base_ring(ode.poly_ring), [p[2] for p in varnames])
    eval_dict =
        Dict(str_to_var(p[1], ode.poly_ring) => str_to_var(p[2], R0) for p in varnames)
    return [eval_at_dict(f, eval_dict) for f in funcs]
end

# -----------------------------------------------------------------------------
