# ------------------------------------------------------------------------------

function total_degree_frac(f::Generic.FracFieldElem{<:MPolyRingElem})
    return sum(map(total_degree, unpack_fraction(f)))
end

function total_degree_frac(f::MPolyRingElem)
    return total_degree(f)
end

# ------------------------------------------------------------------------------

function append_at_index(vec::Vector{T}, idx::Integer, val::T) where {T}
    # NOTE: could also just use insert!
    @assert (idx == 1) || (idx == length(vec) + 1)
    if idx == 1
        vcat(val, vec)
    else
        vcat(vec, val)
    end
end

function cut_at_index(vec::Vector{T}, idx::Integer) where {T}
    @assert (idx == 1) || (idx == length(vec))
    if idx == 1
        vec[2:end]
    else
        vec[1:(end - 1)]
    end
end

# ------------------------------------------------------------------------------

"""
    eval_at_dict(f, d)

Evaluates a polynomial/rational function on a dictionary of type `var => val` and missing values are replaced with zeroes
"""
function eval_at_dict(poly::P, d::Dict{P, <:RingElem}) where {P <: MPolyRingElem}
    R = parent(first(values(d)))
    point = [get(d, v, zero(R)) for v in gens(parent(poly))]
    return evaluate(poly, point)
end

function eval_at_dict(poly::P, d::Dict{P, S}) where {P <: MPolyRingElem, S <: Real}
    R = parent(poly)
    @assert R == parent(first(keys(d)))
    xs = gens(parent(first(keys(d))))
    xs_sym = [get(d, x, 0.0) for x in xs if string(x) in map(string, gens(R))]
    accum = zero(valtype(d))
    for t in terms(poly)
        cf = coeff(t, 1)
        # Nemo.QQ --> Rational{BigInt}
        # NOTE: what about Nemo.GF|ZZ?
        cf_ = BigInt(numerator(cf)) // BigInt(denominator(cf))
        ex = exponent_vector(t, 1)
        accum += cf_ * prod(xs_sym .^ ex)
    end
    return accum
end

function eval_at_dict(
    rational::Generic.FracFieldElem{T},
    d::Dict{T, V},
) where {T <: MPolyRingElem, V}
    f, g = unpack_fraction(rational)
    return eval_at_dict(f, d) / eval_at_dict(g, d)
end

function eval_at_dict(
    rational::Generic.FracFieldElem{<:T},
    d::Dict{T, <:RingElem},
) where {T <: MPolyRingElem}
    f, g = unpack_fraction(rational)
    return eval_at_dict(f, d) * inv(eval_at_dict(g, d))
end

function eval_at_dict(
    rational::Generic.FracFieldElem{<:P},
    d::Dict{<:P, <:Union{<:Generic.FracFieldElem, <:P}},
) where {P <: MPolyRingElem}
    f, g = unpack_fraction(rational)
    return eval_at_dict(f, d) // eval_at_dict(g, d)
end

# ------------------------------------------------------------------------------

function unpack_fraction(f::MPolyRingElem)
    return (f, one(parent(f)))
end

function unpack_fraction(f::Generic.FracFieldElem{<:MPolyRingElem})
    return (numerator(f), denominator(f))
end

# ------------------------------------------------------------------------------

function simplify_frac(numer::P, denom::P) where {P <: MPolyRingElem}
    gcd_sub = gcd(numer, denom)
    sub_numer = divexact(numer, gcd_sub)
    sub_denom = divexact(denom, gcd_sub)
    return sub_numer, sub_denom
end

function is_rational_func_const(f)
    is_constant(numerator(f)) && is_constant(denominator(f))
end

function is_rational_func_normalized(f)
    leading_coefficient(denominator(f)) > 0 && isone(gcd(numerator(f), denominator(f)))
end

"""
    compare_rational_func_by(f, g, by)

Returns 
- `-1` if `f < g`,
- `0` if `f = g`, and 
- `1` if `f > g`.

Functions' numerators and denominators are compared using `by`.
"""
function compare_rational_func_by(f, g, by::By, priority = :numerator) where {By}
    # Specializes on the type of `by`
    numf, denf = unpack_fraction(f)
    numg, deng = unpack_fraction(g)
    keynumf, keydenf = by(numf), by(denf)
    keynumg, keydeng = by(numg), by(deng)
    if priority === :numerator
        keynumf < keynumg && return -1
        keynumf > keynumg && return 1
        keydenf < keydeng && return -1
        keydenf > keydeng && return 1
    elseif priority === :denominator
        keydenf < keydeng && return -1
        keydenf > keydeng && return 1
        keynumf < keynumg && return -1
        keynumf > keynumg && return 1
    elseif priority === :additive
        keydenf + keynumf < keynumg + keydeng && return -1
        keydenf + keynumf > keynumg + keydeng && return 1
    else
        throw(DomainError("Unknown value for keyword argument priority", priority))
    end
    return 0
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
    parent_ring_change(poly, new_ring)

Converts a polynomial to a different polynomial ring
Input
- `poly` - a polynomial to be converted
- `new_ring` - a polynomial ring such that every variable name appearing in poly appears among the generators

Output:
- a polynomial in `new_ring` “equal” to `poly`
"""
function parent_ring_change(
    poly::MPolyRingElem,
    new_ring::MPolyRing;
    matching = :byname,
    shift = 0,
)
    old_ring = parent(poly)
    # Construct a mapping for the variable indices.
    # Zero indicates no image of the old variable in the new ring  
    var_mapping = zeros(Int, max(nvars(old_ring), nvars(new_ring)))
    if matching === :byname
        old_symbols, new_symbols = symbols(old_ring), symbols(new_ring)
        for i in 1:length(old_symbols)
            u = old_symbols[i]
            found = findfirst(v -> (u === v), new_symbols)
            isnothing(found) && continue
            var_mapping[i] = found
        end
    elseif matching === :byindex
        var_mapping[1:(nvars(new_ring) - shift)] .= (1 + shift):nvars(new_ring)
    else
        throw(Base.ArgumentError("Unknown matching type: $matching"))
    end
    # Hoist the compatibility check out of the loop
    for i in 1:nvars(old_ring)
        if degree(poly, i) > 0 && iszero(var_mapping[i])
            throw(
                Base.ArgumentError(
                    """
                    The polynomial $poly contains a variable $(gens(old_ring)[i]) not present in the new ring.
                    New ring variables are $(gens(new_ring)))""",
                ),
            )
        end
    end
    bring = base_ring(new_ring)
    exps = Vector{Vector{Int}}(undef, length(poly))
    coefs = map(c -> bring(c), coefficients(poly))
    @inbounds for i in 1:length(poly)
        evec = exponent_vector(poly, i)
        new_exp = zeros(Int, nvars(new_ring))
        for i in 1:length(evec)
            iszero(var_mapping[i]) && continue
            new_exp[var_mapping[i]] = evec[i]
        end
        exps[i] = new_exp
    end
    return new_ring(coefs, exps)
end

function parent_ring_change(
    f::Generic.FracFieldElem{<:MPolyRingElem},
    new_ring::MPolyRing;
    matching = :byname,
)
    n, d = unpack_fraction(f)
    return parent_ring_change(n, new_ring; matching = matching) //
           parent_ring_change(d, new_ring; matching = matching)
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
    cut_indices = map(v -> findfirst(x -> x == v, xs), variables)
    coeff_indices = setdiff(collect(1:length(xs)), cut_indices)
    coeff_vars = xs[coeff_indices]

    K = base_ring(parent(poly))
    new_ring, _ = Nemo.polynomial_ring(K, map(vv -> var_to_str(vv, xs = xs), coeff_vars))
    FieldType = elem_type(K)

    result = Dict{Vector{Int}, Tuple{Vector{Vector{Int}}, Vector{FieldType}}}()

    @inbounds for i in 1:length(poly)
        coef = coeff(poly, i)
        evec = exponent_vector(poly, i)
        var_slice = [evec[i] for i in cut_indices]
        if !haskey(result, var_slice)
            monom_vect, coef_vect = Vector{Vector{Int}}(), Vector{FieldType}()
            sizehint!(monom_vect, 8)
            sizehint!(coef_vect, 8)
            result[var_slice] = (monom_vect, coef_vect)
        end
        monom_vect, coef_vect = result[var_slice]
        new_monom = Vector{Int}(undef, length(coeff_vars))
        for i in 1:length(new_monom)
            new_monom[i] = evec[coeff_indices[i]]
        end
        push!(monom_vect, new_monom)
        push!(coef_vect, coef)
    end

    return Dict(k => new_ring(v[2], v[1]) for (k, v) in result)
end

# ------------------------------------------------------------------------------

function str_to_var(s::String, ring::MPolyRing)
    ind = findfirst(v -> (string(v) == s), symbols(ring))
    if ind === nothing
        throw(Base.KeyError("Variable $s is not found in ring $ring"))
    end
    return gens(ring)[ind]
end

# ------------------------------------------------------------------------------

function var_to_str(v::MPolyRingElem; xs = gens(parent(v)))
    ind = findfirst(vv -> vv == v, xs)
    return string(symbols(parent(v))[ind])
end

# ------------------------------------------------------------------------------

"""
    gen_tag_name(base; stop_words)
    gen_tag_names(n, base; stop_words)

Generates a string which will not collide with the words in `stop_words`.

## Arguments

- `n`: Generates a sequence of unique strings of length `n`
- `base`: A string or a vector of strings, the base for the generated sequence
- `stop_words`: A vector of strings, stop words
"""
function gen_tag_name(base = "T"; stop_words = Vector{String}())
    return first(gen_tag_names(1, base, stop_words = stop_words))
end

function gen_tag_names(n::Integer, base = "T"; stop_words = Vector{String}())
    sequence = if base isa Vector{String}
        @assert n == length(base)
        base
    else
        repeat([base], n)
    end
    while true
        rand_token = Int(rand(UInt8))
        sequence = map(c -> "$(rand_token)__$c", sequence)
        sequence = map(ic -> "$(ic[2])_$(ic[1])", enumerate(sequence))
        if all(elem -> !(elem in stop_words), sequence)
            break
        end
    end
    return sequence
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

"""
    select_pivots(M::MatElem)

Takes as input a matrix M in the reduced row echelon form and returns
tow lists: of the pivot indices and the non-pivot ones
"""
function select_pivots(M::MatElem)
    @assert is_rref(M)
    j = 1
    (nrows, ncols) = size(M)
    nonpivots = Vector{Int}()
    pivots = Vector{Int}()
    for i in 1:ncols
        pivot = false
        for k in j:nrows
            if !iszero(M[k, i])
                j = k + 1
                pivot = true
            end
        end
        if pivot
            push!(pivots, i)
        else
            push!(nonpivots, i)
        end
    end
    return pivots, nonpivots
end
