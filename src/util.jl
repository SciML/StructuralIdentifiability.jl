# ------------------------------------------------------------------------------

function Nemo.vars(f::Generic.Frac{<:MPolyElem})
    return collect(union(Set(vars(numerator(f))), Set(vars(denominator(f)))))
end

function Nemo.total_degree(f::Generic.Frac{<:MPolyElem})
    return sum(map(total_degree, unpack_fraction(f)))
end

# ------------------------------------------------------------------------------

"""
    eval_at_dict(f, d)

Evaluates a polynomial/rational function on a dictionary of type `var => val` and missing values are replaced with zeroes
"""
function eval_at_dict(poly::P, d::Dict{P, <:RingElem}) where {P <: MPolyElem}
    R = parent(first(values(d)))
    point = [get(d, v, zero(R)) for v in gens(parent(poly))]
    return evaluate(poly, point)
end

function eval_at_dict(poly::P, d::Dict{P, S}) where {P <: MPolyElem, S <: Num}
    R = parent(poly)
    xs = gens(parent(first(keys(d))))
    xs_sym = [d[x] for x in xs if string(x) in map(string, gens(R))]
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

function eval_at_dict(rational::Generic.Frac{T}, d::Dict{T, V}) where {T <: MPolyElem, V}
    f, g = unpack_fraction(rational)
    return eval_at_dict(f, d) // eval_at_dict(g, d)
end

function eval_at_dict(
    rational::Generic.Frac{<:T},
    d::Dict{T, <:RingElem},
) where {T <: MPolyElem}
    f, g = unpack_fraction(rational)
    return eval_at_dict(f, d) * inv(eval_at_dict(g, d))
end

function eval_at_dict(
    rational::Generic.Frac{<:P},
    d::Dict{<:P, <:Union{<:Generic.Frac, <:P}},
) where {P <: MPolyElem}
    f, g = unpack_fraction(rational)
    return eval_at_dict(f, d) // eval_at_dict(g, d)
end

# ------------------------------------------------------------------------------

function unpack_fraction(f::MPolyElem)
    return (f, one(parent(f)))
end

function unpack_fraction(f::Generic.Frac{<:MPolyElem})
    return (numerator(f), denominator(f))
end

# ------------------------------------------------------------------------------

function simplify_frac(numer::P, denom::P) where {P <: MPolyElem}
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
) where {P <: MPolyElem}
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
- a polynomial in `new_ring` “equal” to `poly`
"""
function parent_ring_change(poly::MPolyElem, new_ring::MPolyRing; matching = :byname)
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
        var_mapping[1:nvars(new_ring)] .= 1:nvars(new_ring)
    else
        throw(Base.ArgumentError("Unknown matching type: $matching"))
    end
    # Hoist the compatibility check out of the loop
    for i in 1:nvars(old_ring)
        if !iszero(degree(poly, i))
            if iszero(var_mapping[i])
                throw(
                    Base.ArgumentError(
                        "The polynomial contains a variable $(gens(old_ring)[i]) not present in the new ring $poly",
                    ),
                )
            end
        end
    end
    bring = base_ring(new_ring)
    exps = Vector{Vector{Int}}(undef, length(poly))
    coefs = map(c -> bring(c), coefficients(poly))
    @inbounds for i in 1:length(poly)
        evec = exponent_vector(poly, i)
        new_exp = zeros(Int, nvars(new_ring))
        for i in 1:length(evec)
            new_exp[var_mapping[i]] = evec[i]
        end
        exps[i] = new_exp
    end
    return new_ring(coefs, exps)
end

function parent_ring_change(
    f::Generic.Frac{<:MPolyElem},
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
function uncertain_factorization(f::MPolyElem{fmpq})
    vars_f = vars(f)
    if isempty(vars_f)
        return Array{Tuple{typeof(f), Bool}, 1}()
    end
    main_var = vars_f[end]
    d = Nemo.degree(f, main_var)
    lc_f = coeff(f, [main_var], [d])
    gcd_coef = lc_f
    for i in (d - 1):-1:0
        gcd_coef = gcd(gcd_coef, coeff(f, [main_var], [i]))
    end
    f = divexact(f, gcd_coef)
    lc_f = coeff(f, [main_var], [d])

    is_irr = undef
    while true
        plugin = rand(5:10, length(vars_f) - 1)
        if evaluate(lc_f, vars_f[1:(end - 1)], plugin) != 0
            f_sub = evaluate(f, vars_f[1:(end - 1)], plugin)
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
    runtime = @elapsed prelim_factors = uncertain_factorization(poly)
    _runtime_logger[:id_uncertain_factorization] += runtime
    cert_factors = map(pair -> pair[1], filter(f -> f[2], prelim_factors))
    uncert_factors = map(pair -> pair[1], filter(f -> !f[2], prelim_factors))
    for p in uncert_factors
        for f in Nemo.factor(p)
            push!(cert_factors, f[1])
        end
    end
    push!(_runtime_logger[:id_certain_factors], cert_factors)
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
- `variables` - a list of variables from the generators of the ring of p
Output:
- dictionary with keys being tuples of length `lenght(variables)` and values being polynomials in the variables other than those which are the coefficients at the corresponding monomials (in a smaller polynomial ring)
"""
function extract_coefficients(poly::P, variables::Array{P, 1}) where {P <: MPolyElem}
    xs = gens(parent(poly))
    @assert all(in(xs), variables)
    cut_indices = map(v -> findfirst(x -> x == v, xs), variables)
    coeff_indices = setdiff(collect(1:length(xs)), cut_indices)
    coeff_vars = xs[coeff_indices]

    K = base_ring(parent(poly))
    new_ring, _ = Nemo.PolynomialRing(K, map(vv -> var_to_str(vv, xs = xs), coeff_vars))
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
    if ind == nothing
        throw(Base.KeyError("Variable $s is not found in ring $ring"))
    end
    return gens(ring)[ind]
end

# ------------------------------------------------------------------------------

function var_to_str(v::MPolyElem; xs = gens(parent(v)))
    ind = findfirst(vv -> vv == v, xs)
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
    e = Symbolics.value(e)
    return eval_at_nemo(e, vals)
end

function eval_at_nemo(e::SymbolicUtils.BasicSymbolic, vals::Dict)
    if Symbolics.istree(e)
        # checking if it is a function of the form x(t), a bit dirty
        if length(Symbolics.arguments(e)) == 1 && "$(first(Symbolics.arguments(e)))" == "t"
            return vals[e]
        end
        # checking if this is a vector entry like x(t)[1]
        if Symbolics.operation(e) == getindex
            return vals[e]
        end
        # otherwise, this is a term
        args = map(a -> eval_at_nemo(a, vals), Symbolics.arguments(e))
        if Symbolics.operation(e) in (+, -, *)
            return Symbolics.operation(e)(args...)
        elseif isequal(Symbolics.operation(e), /)
            return //(args...)
        elseif isequal(Symbolics.operation(e), ^)
            if args[2] >= 0
                return args[1]^args[2]
            end
            return 1 // args[1]^(-args[2])
        end
        throw(Base.ArgumentError("Function $(Symbolics.operation(e)) is not supported"))
    elseif e isa Symbolics.Symbolic
        return get(vals, e, e)
    end
end

function eval_at_nemo(e::Union{Integer, Rational}, vals::Dict)
    return e
end

function eval_at_nemo(e::Union{Float16, Float32, Float64}, vals::Dict)
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
If yes, returns a pair (a, number), otherwise nothing
"""
function decompose_derivative(varname::String, prefixes::Array{String})
    for pr in prefixes
        if startswith(varname, pr) && length(varname) > length(pr) + 1
            if varname[length(pr) + 1] == '_' &&
               all(map(isdigit, collect(varname[(length(pr) + 2):end])))
                return (pr, parse(Int, varname[(length(pr) + 2):end]))
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

function get_measured_quantities(ode::ModelingToolkit.ODESystem)
    if any(ModelingToolkit.isoutput(eq.lhs) for eq in ModelingToolkit.equations(ode))
        @info "Measured quantities are not provided, trying to find the outputs in input ODE."
        return filter(
            eq -> (ModelingToolkit.isoutput(eq.lhs)),
            ModelingToolkit.equations(ode),
        )
    else
        throw(
            error(
                "Measured quantities (output functions) were not provided and no outputs were found.",
            ),
        )
    end
end
