import GroebnerBasis

#------------------------------------------------------------------------------

function eval_at_dict(poly, d)
    """
    Evaluates a polynomial on a dict var => val
    missing values are replaced with zeroes
    """
    point = [get(d, v, base_ring(parent(poly))(0)) for v in gens(parent(poly))]
    return evaluate(poly, point)
end

#------------------------------------------------------------------------------

function unpack_fraction(f)
    """
    Maps polynomial/rational function to a pair denominator/numerator
    """
    if applicable(numerator, f)
        return (numerator(f), denominator(f))
    end
    return (f, parent(f)(1))
end

#------------------------------------------------------------------------------

function simplify_frac(numer, denom)
    gcd_sub = gcd(numer, denom)
    sub_numer = divexact(numer, gcd_sub)
    sub_denom = divexact(denom, gcd_sub)
    return sub_numer, sub_denom
end

#------------------------------------------------------------------------------

function make_substitution(f, var_sub, val_numer, val_denom)
    """
    Substitute a variable in a polynomial with an expression
    Input:
        - f, the polynomial
        - var_sub, the variable to be substituted
        - var_numer, numerator of the substitution expression
        - var_denom, denominator of the substitution expression
    Output:
        polynomial, result of substitution
    """
    d = degree(f, var_sub)
    f_subs = sum(coeff(f, [var_sub], [i]) * (val_numer ^ i) * (val_denom ^ (d - i)) for i in 0:d)
    return f_subs
end

#------------------------------------------------------------------------------

function evaluate_frac(f, var_sub, val_numer, val_denom)
    f_numer = make_substitution(numerator(f), var_sub, val_numer, val_denom) // val_denom^(degree(numerator(f), var_sub))
    f_denom = make_substitution(denominator(f), var_sub, val_numer, val_denom) // val_denom^(degree(denominator(f), var_sub))
    return f_numer // f_denom
end

#------------------------------------------------------------------------------

function parent_ring_change(poly, new_ring)
    """
    Converts a polynomial to a different polynomial ring
    Input
      - poly - a polynomial to be converted
      - new_ring - a polynomial ring such that every variable name
          appearing in poly appears among the generators
    Output: a polynomial in new_ring "equal" to poly
    """
    old_ring = parent(poly)
    # construct a mapping for the variable indices
    var_mapping = []
    for u in gens(old_ring)
        push!(
            var_mapping,
            findfirst(v -> (string(u) == string(v)), gens(new_ring))
        )
    end
    builder = MPolyBuildCtx(new_ring)
    for term in zip(exponent_vectors(poly), coeffs(poly))
        exp, coef = term
        new_exp = [0 for _ in gens(new_ring)]
        for i in 1:length(exp)
            if exp[i] != 0
                if var_mapping[i] == nothing
                    throw(Base.ArgumentError("The polynomial contains a variable not present in the new ring $poly"))
                else
                    new_exp[var_mapping[i]] = exp[i]
                end
            end
        end
        push_term!(builder, new_ring.base_ring(coef), new_exp)
    end
    return finish(builder)
end

#------------------------------------------------------------------------------

function check_factors(f)
    """
    !!! gens(parent(f))[end] must be y_order (highest order of y)!
    If return (true, 1), f must be irreducible.
    If return (true, gcd_coef), gcd_coef is an extra factor of f, and f//gcd_coef is irreducible.
    If return (false, 1), f may or may not be irredicible.
    If return (false, gcd_coef), gcd_coef is an extra factor of f, f//gcd_coef may or may not be irredicible.
    """
    vars_f = vars(f)
    var_test = vars_f[end]
    d = degree(f, var_test)
    lc = coeff(f, [var_test], [d])
    is_irr = true
    while true
        plugin = rand(5:10, length(vars_f) - 1)
        if evaluate(lc, vars_f[1 : end - 1], plugin) != 0    
            f_sub = evaluate(f, vars_f[1 : end - 1], plugin)
            uni_ring, var_uni = PolynomialRing(base_ring(f), "v")
            f_uni = to_univariate(uni_ring, f_sub)
            is_irr = isirreducible(f_uni)
            break
        end
    end
    gcd_coef = lc
    for i in  (d - 1):-1:0
        gcd_coef = gcd(gcd_coef, coeff(f, [var_test], [i]))
    end
    return (is_irr, gcd_coef)
end

#------------------------------------------------------------------------------

function dict_to_poly(dict_monom, poly_ring)
    builder = MPolyBuildCtx(poly_ring)
    for (monom, coef) in pairs(dict_monom)
        push_term!(builder, poly_ring.base_ring(coef), monom)
    end
    return finish(builder)
end

#------------------------------------------------------------------------------

function extract_coefficients(poly, variables)
    """
    Intput:
        poly - multivariate polynomial
        variables - a list of variables from the generators of the ring of p
    Output:
        dictionary with keys being tuples of length len(variables) and values being 
        polynomials in the variables other than variables which are the coefficients
        at the corresponding monomials (in a smaller polynomial ring)
    """
    var_to_ind = Dict([(v, findfirst(e -> (e == v), gens(parent(poly)))) for v in variables])
    indices = [var_to_ind[v] for v in variables]

    coeff_vars = filter(v -> !(string(v) in map(string, variables)), gens(parent(poly)))
    new_ring, new_vars = PolynomialRing(base_ring(parent(poly)), map(string, coeff_vars))
    coeff_var_to_ind = Dict([(v, findfirst(e -> (e == v), gens(parent(poly)))) for v in coeff_vars])

    result = Dict()

    for (monom, coef) in zip(exponent_vectors(poly), coeffs(poly))
        var_slice = [monom[i] for i in indices]
        if !haskey(result, var_slice)
            result[var_slice] = Dict()
        end
        new_monom = [0 for _ in 1:length(coeff_vars)]
        for i in 1:length(new_monom)
            new_monom[i] = monom[coeff_var_to_ind[coeff_vars[i]]]
        end
        result[var_slice][new_monom] = coef
    end

    return Dict([(k, dict_to_poly(v, new_ring)) for (k, v) in pairs(result)])
end

#------------------------------------------------------------------------------

function check_injectivity(polys; method="Singular")
    """
    Checks a generic injectivity of the *projective* map defined by polys
    Inputs:
        - polys - a list of polynomials
    Output: dictionary from variables of the parent ring to booleans as follows:
        True: every generic fiber of the map defined by polys has a single value of the variabe
        False: otherwise
    """
    @debug "Constructing equations"
    flush(stdout)
    ring = parent(polys[1])
    point = map(v -> rand(5:100), gens(ring))
    pivot = sort(polys, by = (p -> total_degree(ring(p))))[1]
    @debug "Pivot polynomial is $pivot"
    flush(stdout)
    eqs = []
    for p = polys
        push!(eqs, p * evaluate(ring(pivot), point) - evaluate(ring(p), point) * pivot)
    end
    ring_sing, vars_sing = Singular.PolynomialRing(Singular.QQ, vcat(map(string, gens(ring)), "sat_aux"); ordering=:degrevlex)
    eqs_sing = map(p -> parent_ring_change(p, ring_sing), eqs)
    push!(
        eqs_sing,
        parent_ring_change(pivot, ring_sing) * vars_sing[end] - 1
    )
    @debug "Computing Groebner basis ($(length(eqs_sing)) equations)"
    flush(stdout)
    if method == "Singular"
        gb = Singular.std(Singular.Ideal(ring_sing, eqs_sing))
    elseif method == "GroebnerBasis"
        gb = GroebnerBasis.f4(Singular.Ideal(ring_sing, eqs_sing))
    else
        throw(Base.ArgumentError("Unknown method $method"))
    end

    @debug "Producing the result"
    flush(stdout)
    result = Dict()
    for v in gens(ring)
        v_sing = parent_ring_change(v, ring_sing)
        result[v] = (total_degree(Singular.reduce(v_sing, gb)) == 0)
    end
    return result
end

#------------------------------------------------------------------------------

function check_identifiability(io_equation, parameters; method="Singular")
    """
    For the io_equation and the list of all parameter variables, returns a dictionary
    var => whether_globally_identifiable
    method can be "Singular" or "GroebnerBasis" yielding using Singular.jl or GroebnerBasis.jl
    """
    @debug "Extracting coefficients"
    flush(stdout)
    nonparameters = filter(v -> !(string(v) in map(string, parameters)), gens(parent(io_equation)))
    coeffs = extract_coefficients(io_equation, nonparameters)

    return check_injectivity(collect(values(coeffs)); method=method)
end

#------------------------------------------------------------------------------

function str_to_var(s, ring)
    return gens(ring)[findfirst(v -> (string(v) == s), gens(ring))]
end

#------------------------------------------------------------------------------
