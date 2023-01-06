# ------------------------------------------------------------------------------

function Nemo.vars(f::Generic.Frac{<: MPolyElem})
    return collect(union(Set(vars(numerator(f))), Set(vars(denominator(f)))))
end

# ------------------------------------------------------------------------------

"""
    eval_at_dict(f, d)

Evaluates a polynomial/rational function on a dictionary of type `var => val` and missing values are replaced with zeroes
"""
function eval_at_dict(poly::P, d::Dict{P,<: RingElem}) where P <: MPolyElem
    R = parent(first(values(d)))
    point = [get(d, v, zero(R)) for v in gens(parent(poly))]
    return evaluate(poly, point)
end

function eval_at_dict(rational::Generic.Frac{<: P}, d::Dict{P,<: RingElem}) where P <: MPolyElem
    f, g = unpack_fraction(rational)
    return eval_at_dict(f, d) * inv(eval_at_dict(g, d))
end

function eval_at_dict(rational::Generic.Frac{<: P}, d::Dict{P,<: Generic.Frac}) where P <: MPolyElem
    f, g = unpack_fraction(rational)
    return eval_at_dict(f, d) // eval_at_dict(g, d)
end

function eval_at_dict(rational::Generic.Frac{<: P}, d::Dict{P,<: MPolyElem}) where P <: MPolyElem
    f, g = unpack_fraction(rational)
    return eval_at_dict(f, d) // eval_at_dict(g, d)
end

# ------------------------------------------------------------------------------

function unpack_fraction(f::MPolyElem)
    return (f, one(parent(f)))
end

function unpack_fraction(f::Generic.Frac{<: MPolyElem})
    return (numerator(f), denominator(f))
end

# ------------------------------------------------------------------------------

function simplify_frac(numer::P, denom::P) where P <: MPolyElem
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
function make_substitution(f::P, var_sub::P, val_numer::P, val_denom::P) where P <: MPolyElem
    d = Nemo.degree(f, var_sub)

    result = 0
    @debug "Substitution in a polynomial of degree $d"
    flush(stdout)
    for i in 0:d
        @debug "\t Degree $i"
        flush(stdout)
        result += coeff(f, [var_sub], [i]) * (val_numer^i) * (val_denom^(d - i))
        @debug "\t Intermediate result of size $(length(result))"
    end
    return result
end

# ------------------------------------------------------------------------------

"""
    parent_ring_change(poly, new_ring)

Converts a polynomial to a different polynomial ring
Input
- `poly` - a polynomial to be converted
- `new_ring` - a polynomial ring such that every variable name appearing in poly appears among the generators

Output:
- a polynomial in `new_ring` "equal" to `poly`
"""
function parent_ring_change(poly::MPolyElem, new_ring::MPolyRing; matching=:byname)
    old_ring = parent(poly)
    # construct a mapping for the variable indices
    var_mapping = Array{Any,1}()

    if matching == :byname
        for u in symbols(old_ring)
            push!(
                var_mapping,
                findfirst(v -> (string(u) == string(v)), symbols(new_ring))
            )
        end
    elseif matching == :byindex
        append!(var_mapping, 1:length(symbols(new_ring)))
        if length(symbols(new_ring)) < length(symbols(old_ring))
            append!(var_mapping, Array{Any,1}[nothing, length(symbols(old_ring)) - length(symbols(new_ring))])
        end
    else
        throw(Base.ArgumentError("Unknown matching type: $matching"))
    end
    builder = MPolyBuildCtx(new_ring)
    for term in zip(exponent_vectors(poly), coefficients(poly))
        exp, coef = term
        new_exp = [0 for _ in gens(new_ring)]
        for i in 1:length(exp)
            if exp[i] != 0
                if var_mapping[i] == nothing
                    throw(Base.ArgumentError("The polynomial contains a variable $(gens(old_ring)[i]) not present in the new ring $poly"))
                else
                    new_exp[var_mapping[i]] = exp[i]
                end
            end
        end
        push_term!(builder, new_ring.base_ring(coef), new_exp)
    end
    return finish(builder)
end

function parent_ring_change(f::Generic.Frac{<: MPolyElem}, new_ring::MPolyRing; matching=:byname)
    n, d = unpack_fraction(f)
    return parent_ring_change(n, new_ring; matching=matching) // parent_ring_change(d, new_ring; matching=matching)
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
function uncertain_factorization(f::MPolyElem{fmpq})
    vars_f = vars(f)
    if isempty(vars_f)
        return Array{Tuple{typeof(f),Bool},1}()
    end
    main_var = vars_f[end]
    d = Nemo.degree(f, main_var)
    lc_f = coeff(f, [main_var], [d])
    gcd_coef = lc_f
    for i in  (d - 1):-1:0
        gcd_coef = gcd(gcd_coef, coeff(f, [main_var], [i]))
    end
    f = divexact(f, gcd_coef)
    lc_f = coeff(f, [main_var], [d])

    is_irr = undef
    while true
        plugin = rand(5:10, length(vars_f) - 1)
        if evaluate(lc_f, vars_f[1:end - 1], plugin) != 0
            f_sub = evaluate(f, vars_f[1:end - 1], plugin)
            uni_ring, var_uni = Nemo.PolynomialRing(base_ring(f), string(main_var))
            f_uni = to_univariate(uni_ring, f_sub)
            # if !issquarefree(f_uni)
            if !isone(gcd(f, derivative(f, main_var)))
                @debug "GCD is $(gcd(f, derivative(f, main_var)))"
                f = divexact(f, gcd(f, derivative(f, main_var)))
            end
            # end
            is_irr = Nemo.isirreducible(f_uni)
            break
        end
    end

    coeff_factors = uncertain_factorization(gcd_coef)
    push!(coeff_factors, (f, is_irr))
end

# ------------------------------------------------------------------------------

function fast_factor(poly::MPolyElem{fmpq})
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

function dict_to_poly(dict_monom::Dict{Array{Int,1},<: RingElem}, poly_ring::MPolyRing)
    builder = MPolyBuildCtx(poly_ring)
    for (monom, coef) in pairs(dict_monom)
        push_term!(builder, poly_ring.base_ring(coef), monom)
    end
    return finish(builder)
end

# ------------------------------------------------------------------------------

"""
    extract_coefficients(poly, variables)

Intput:
- `poly` - multivariate polynomial
- `variables` - a list of variables from the generators of the ring of p
Output:
- dictionary with keys being tuples of length `lenght(variables)` and values being polynomials in the variables other than those which are the coefficients at the corresponding monomials (in a smaller polynomial ring)
"""
function extract_coefficients(poly::P, variables::Array{P,1}) where P <: MPolyElem
    var_to_ind = Dict((v, findfirst(e -> (e == v), gens(parent(poly)))) for v in variables)
    indices = [var_to_ind[v] for v in variables]

    coeff_vars = filter(v -> !(var_to_str(v) in map(var_to_str, variables)), gens(parent(poly)))
    new_ring, new_vars = Nemo.PolynomialRing(base_ring(parent(poly)), map(var_to_str, coeff_vars))
    coeff_var_to_ind = Dict((v, findfirst(e -> (e == v), gens(parent(poly)))) for v in coeff_vars)
    FieldType = typeof(one(base_ring(new_ring)))

    result = Dict{Array{Int,1},Dict{Array{Int,1},FieldType}}()

    for (monom, coef) in zip(exponent_vectors(poly), coefficients(poly))
        var_slice = [monom[i] for i in indices]
        if !haskey(result, var_slice)
            result[var_slice] = Dict{Array{Int,1},FieldType}()
        end
        new_monom = [0 for _ in 1:length(coeff_vars)]
        for i in 1:length(new_monom)
            new_monom[i] = monom[coeff_var_to_ind[coeff_vars[i]]]
        end
        result[var_slice][new_monom] = coef
    end

    return Dict(k => dict_to_poly(v, new_ring) for (k, v) in result)
end

# ------------------------------------------------------------------------------

function str_to_var(s::String, ring::MPolyRing)
    ind = findfirst(v -> (string(v) == s), symbols(ring))
    if ind == nothing
        throw(Base.KeyError("Variable $s is not found in ring $ring"))
    end
    return gens(ring)[ind]
end

# ------------------------------------------------------------------------------

function var_to_str(v::MPolyElem)
    ind = findfirst(vv -> vv == v, gens(parent(v)))
    return string(symbols(parent(v))[ind])
end

# ------------------------------------------------------------------------------

"""
    switch_ring(v, ring)

For a variable `v`, returns a variable in `ring` with the same name
"""
function switch_ring(v::MPolyElem, ring::MPolyRing)
    ind = findfirst(vv -> vv == v, gens(parent(v)))
    return str_to_var(string(symbols(parent(v))[ind]), ring)
end

# ------------------------------------------------------------------------------

function eval_at_nemo(e::Num, vals::Dict)
    e = ModelingToolkit.Symbolics.value(e)
    return eval_at_nemo(e, vals)
end

function eval_at_nemo(e, vals::Dict)
    if ModelingToolkit.Symbolics.istree(e)
        args = map(a -> eval_at_nemo(a, vals), ModelingToolkit.Symbolics.arguments(e))
        if ModelingToolkit.Symbolics.operation(e) in [+, -, *]
            return ModelingToolkit.Symbolics.operation(e)(args...)
        elseif isequal(ModelingToolkit.Symbolics.operation(e), /)
            return //(args...)
        end
        if ModelingToolkit.Symbolics.operation(e) === ^
            if args[2] >= 0
                return args[1]^args[2]
            end
            return 1 // args[1]^(-args[2])
        end
        throw(Base.ArgumentError("Function $(ModelingToolkit.Symbolics.operation(e)) is not supported"))
    end
end

function eval_at_nemo(e::Union{ModelingToolkit.Symbolics.Sym,ModelingToolkit.Symbolics.Term}, vals::Dict)
    if typeof(e) <: ModelingToolkit.Symbolics.Term{Real,Nothing}
        throw(Base.ArgumentError("Function $(ModelingToolkit.Symbolics.operation(e)) is not supported"))
    end
    return get(vals, e, e)
end

function eval_at_nemo(e::Union{Integer,Rational}, vals::Dict)
    return e
end

function eval_at_nemo(e::Union{Float16,Float32,Float64}, vals::Dict)
    if isequal(e % 1, 0)
        out = Int(e)
    else
        out = rationalize(e)
    end
    @warn "Floating points are not allowed, value $e will be converted to $(out)."
    return out
end

# -----------------------------------------------------------------------------

"""
    decompose_derivative(varname, prefixes)

Determines if it is possible to represent the `varname` as `a_number` where `a` is an element of `prefixes`
"""
function decompose_derivative(varname::String, prefixes::Array{String})
    for pr in prefixes
        if startswith(varname, pr) && length(varname) > length(pr) + 1
            if varname[length(pr) + 1] == '_' && all(map(isdigit, collect(varname[length(pr) + 2:end])))
                return (pr, parse(Int, varname[length(pr) + 2:end])) 
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

function difforder(diffpoly::MPolyElem, prefix::String)
    orders = [-1]
    for v in vars(diffpoly)
        d = decompose_derivative(var_to_str(v), [prefix])
        if d != nothing
            push!(orders, d[2])
        end
    end
    return max(orders...)
end
