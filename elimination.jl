using Dates
using Logging
using Oscar

include("util.jl")

#------------------------------------------------------------------------------

function det_minor_expansion_inner(m, discarded, cache)
    n = size(m, 1);
    if length(discarded[1]) == n
        return 1;
    end
    if discarded in keys(cache)
        return cache[discarded];
    end
    result = 0;
    row = minimum(setdiff(Set(1:n), Set(discarded[1])));
    dis_rows = Tuple(sort([[i for i in discarded[1]]; row]));
    sign = 1;
    for col = 1:n
        if !(col in discarded[2])
            dis_cols = Tuple(sort([[i for i in discarded[2]]; col]));
            result += sign * m[row, col] * det_minor_expansion_inner(m, (dis_rows, dis_cols), cache);
            sign = -1 * sign;
        end
    end
    if length(discarded[1]) > 1
        cache[discarded] = result;
    end
    if length(discarded[1]) < 3
        @debug "Discarded: $discarded; $(Dates.now())"
    end
    return result;
end

function det_minor_expansion(m)
    cache = Dict();
    return det_minor_expansion_inner(m, (Tuple{}(), Tuple{}()), cache);
end

#------------------------------------------------------------------------------

function Bezout_matrix(f, g, var_elim)
    """
    Compute the Bezout matrix of two polynomials f, g with respect to var_elim
    Inputs:
        - f::AbstractAlgebra.MPolyElem, first polynomial
        - g::AbstractAlgebra.MPolyElem, second polynomial
        - var_elim::AbstractAlgebra.MPolyElem, variable, of which f and g are considered as polynomials
    Output:
        - M::MatrixElem, The Bezout matrix
    """
    parent_ring = parent(f)
    deg_f = degree(f, var_elim)
    deg_g = degree(g, var_elim)
    n = max(deg_f, deg_g)
    coeffs_f = [coeff(f, [var_elim], [i]) for i in 0:n]
    coeffs_g = [coeff(g, [var_elim], [i]) for i in 0:n]
    GL = MatrixSpace(parent_ring, n, n)
    M = zero(GL)
    for i in 1:n
        for j in 1:n
            M[i, j] = sum([coeffs_f[j + k + 1] * coeffs_g[i - k] - coeffs_g[j + k + 1] * coeffs_f[i - k] for k in 0:min(i - 1, n - j)])
        end
    end
    return M
end

#------------------------------------------------------------------------------

function simplify_matrix(M)
    """
    Eliminate GCD of entries of every row and column
    Input:
        - M::MatrixElem, matrix to be simplified
    Output:
        - M::MatrixElem, Simplified matrix
        - extra_factors::Vector{AbstractAlgebra.MPolyElem}, array of GCDs eliminated from M.
    """

    function _simplify_range(coords)
        """
        An auxiliary function taking a list of coordinates of cells
        and dividing them by their gcd.
        Returns the gcd
        """
        gcd_temp = M[coords[1]...]
        for c in coords[2:end]
            gcd_temp = gcd(gcd_temp, M[c...])
        end
        if gcd_temp != 1
            for c in coords
                M[c...] = divexact(M[c...], gcd_temp)
            end
        end
        return gcd_temp
    end

    extra_factors = []
    rows_cols = []
    # adding all rows
    for i in 1:nrows(M)
        push!(rows_cols, [(i, j) for j in 1:ncols(M)])
    end
    # adding all columns
    for j in 1:ncols(M)
        push!(rows_cols, [(i, j) for i in 1:nrows(M)])
    end
    # eliminating factors
    for r in rows_cols
        gcd_temp = _simplify_range(r)
        if gcd_temp != 1
            push!(extra_factors, gcd_temp)
        end
    end

    return M, extra_factors
end

#------------------------------------------------------------------------------

# Definition: for an irreducible polynomial P in n variables, we will call 
# an iterator Julia object *generic point generator* if it generates 
# an infinite random sequence of points in C^n such that
#    1. P vanishes at each of the points
#    2. for every other irreducible polynomial, the probability of
#    vanishing at each of these points is zero
#
# Each returned point is dictionary from the variables to the values

#------------------------------------------------------------------------------

mutable struct RationalVarietyPointGenerator
    ring
    equations
    parametric_vars
    ind_to_nonparam
    linear_system_A
    linear_system_b
    cached_points
    function RationalVarietyPointGenerator(equations, parametric_vars)
        ring = parent(equations[1])
        nonparametric_vars = filter(v -> !(v in parametric_vars), gens(ring))
        ind_to_nonparam = Dict(i => nonparametric_vars[i] for i in 1:length(nonparametric_vars))
        codim = length(nonparametric_vars)
        S = MatrixSpace(ring, codim, codim)
        Sv = MatrixSpace(ring, codim, 1)

        linear_system_A = zero(S)
        linear_system_b = zero(Sv)

        for i in 1:codim
            reminder = equations[i]
            for j in 1:codim
                linear_system_A[i, j] = derivative(equations[i], nonparametric_vars[j])
                reminder = reminder - linear_system_A[i, j] * nonparametric_vars[j]
            end
            linear_system_b[i, 1] = -reminder
        end

        return new(ring, equations, parametric_vars, ind_to_nonparam, linear_system_A, linear_system_b, [])
    end
end

#------------------------------------------------------------------------------

function Base.iterate(p::RationalVarietyPointGenerator, i::Int=1)
    if i > length(p.cached_points)
        @debug "Generating new point on the variety"
        sample_max = i * 50
        result = undef
        while true
            result = Dict(v => base_ring(p.ring)(rand(1:sample_max)) for v in p.parametric_vars)
            A = map(poly -> eval_at_dict(poly, result), p.linear_system_A)
            b = map(poly -> eval_at_dict(poly, result), p.linear_system_b)
            x = undef
            try
                x = solve(A, b)
            catch e
                continue
            end
            for i in 1:length(x)
                result[p.ind_to_nonparam[i]] = x[i]
            end
            break
        end
        push!(p.cached_points, result)
    end
    return (p.cached_points[i], i + 1)
end

#------------------------------------------------------------------------------

function choose(polys, generic_point_generator)
    """
    Input:
        - array_f, an array of distinct irreducible polynomials in the same ring
        - generic_point_generator, a generic point generator as described above for one of polys
    Output:
        - the polynomial that vanishes at the generic_point_generator
    """
    vars = gens(parent(polys[1]))
    for p in generic_point_generator
        if length(polys) <= 1
            break
        end
        point = [p[v] for v in vars]
        filter!(e -> (evaluate(e, point) == 0), polys)
    end
    return polys[1]
end

#------------------------------------------------------------------------------

function eliminate_var(f, g, var_elim, generic_point_generator)
    """
    Eliminate variable from a pair of polynomials
    Input:
        - f and g, polynomials
        - var_elim, variable to be eliminated
        - generic_point_generator, a generic point generator object for the factor
          of the resultant of f and g of interest
    Output:
        polynomial, the desired factor of the resultant of f and g
    """
    #Step 1: Possible simplification for (f,g)
    #Linear comb
    while f != 0 && g != 0
        if degree(f, var_elim) > degree(g, var_elim)
            f, g = g, f
        end
        lf = coeff(f, [var_elim], [degree(f, var_elim)])
        lg = coeff(g, [var_elim], [degree(g, var_elim)])
        (flag, q) = divides(lg, lf)
        if flag
            @debug "\t Decreasing degree with linear combination $(Dates.now())"
            flush(stdout)
            g = g - q * f * var_elim ^ (degree(g, var_elim) - degree(f, var_elim))
        elseif (degree(g, var_elim) == degree(f, var_elim))
            (flag, q) = divides(lf, lg)
            if flag
                @debug "\t Decreasing degree with linear combination $(Dates.now())"
                flush(stdout)
                f = f - q * g
            else
                break
            end
        else
            break
        end
    end
    
    if f == 0 || g == 0
        return f + g
    end

    #Case (f(v^d), g(v^d)):
    list_deg = []
    for d in 0:degree(f, var_elim)
        if coeff(f, [var_elim], [d]) != 0
            push!(list_deg, d)
        end
    end
    for d in 0:degree(g, var_elim)
        if coeff(g, [var_elim], [d]) != 0
            push!(list_deg, d) 
        end
    end
    gcd_deg = list_deg[1]
    for ele in list_deg[2:end]
        gcd_deg = gcd(gcd_deg, ele)
    end
    if gcd_deg > 1
        f = sum([coeff(f, [var_elim], [gcd_deg * i]) * (var_elim ^ i) for i in 0:(degree(f, var_elim) ÷ gcd_deg)])
        g = sum([coeff(g, [var_elim], [gcd_deg * i]) * (var_elim ^ i) for i in 0:(degree(g, var_elim) ÷ gcd_deg)])
    end

    #Step 2: Initialization
    extra_factors = []

    #Step 3: Compute resultant
    if degree(f, var_elim) == 0
        R = f
    elseif degree(f, var_elim) == 1
        coef_1 = coeff(f, [var_elim], [1])
        coef_0 = coeff(f, [var_elim], [0])
        R = make_substitution(g, var_elim, -coef_0, coef_1)
    else
        @debug "Calculating Bezout Matrix $(Dates.now())"
        flush(stdout)
        M = Bezout_matrix(f, g, var_elim)
        if generic_point_generator != nothing
            @debug "Simplifying Bezout Matrix $(Dates.now())"
            flush(stdout)
            M_simp, extra_factors = simplify_matrix(M)
        else
            M_simp = M
        end
        M_size = zero(MatrixSpace(ZZ, ncols(M_simp), ncols(M_simp)))
        for i in 1:ncols(M_simp)
            for j in 1:ncols(M_simp)
                M_size[i,j] = length(M_simp[i,j])
            end
        end
        @debug "\t Bezout matrix size: \n $M_size"
        @debug "\t Computing determinant $(Dates.now())"
        flush(stdout)
        R = det_minor_expansion(M_simp)
    end
    #Step 4: Eliminate extra factors
    if generic_point_generator != nothing
        #Preliminary factorization
        #TODO: Theoretical justification
        is_irr, gcd_coef = check_factors(R)
        if is_irr
            @debug "\t Using GCD to eliminate extra factors $(Dates.now())" 
            flush(stdout)
            res_pre = divexact(R, gcd_coef)
            push!(extra_factors, gcd_coef)
            push!(extra_factors, res_pre)
            res = choose(extra_factors, generic_point_generator)
            if res == res_pre
                if gcd_coef != 1
                    @debug "\t \t Size of extra factor: $(length(gcd_coef)); $(Dates.now())"
                    flush(stdout)
                end
                return res
            end
        end
        @debug "\t Preliminary factorization failed, using Singular; $(Dates.now()) "
        flush(stdout)
        var_names = [string(v) for v in gens(parent(f))]
        ring_sing, var_sing = Singular.PolynomialRing(Singular.QQ, var_names) 
        poly_sing = mpoly_conversion(R, ring_sing)        
        R_factors = Singular.factor(poly_sing)
        @debug "\t Size and multiplicity of factors; $(Dates.now()) "
        flush(stdout)
        for fac in R_factors
            fac_Nemo = mpoly_conversion(fac[1], parent(f))
            @debug "Size $(length(fac_Nemo)) -- $(fac[2]) times; $(Dates.now())"
            flush(stdout)
            push!(extra_factors, fac_Nemo)
        end
        res = choose(extra_factors, generic_point_generator)
    else
        res = R
    end
    return res
end
