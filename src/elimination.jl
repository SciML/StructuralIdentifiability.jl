# ------------------------------------------------------------------------------

PairIntTuples = Tuple{Tuple{Vararg{Int}},Tuple{Vararg{Int}}}

function det_minor_expansion_inner(
        m::MatElem{<: T}, 
        discarded::PairIntTuples, 
        cache::Dict{PairIntTuples,T}
    ) where T <: RingElem
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
        flush(stdout)
    end
    return result;
end

function det_minor_expansion(m::MatElem{T}) where T <: RingElem
    cache = Dict{PairIntTuples,T}();
    return det_minor_expansion_inner(m, (Tuple{}(), Tuple{}()), cache);
end

# ------------------------------------------------------------------------------

"""
    Bezout_matrix(f, g, var_elim)

Compute the Bezout matrix of two polynomials f, g with respect to var_elim
Inputs:
    - f - first polynomial
    - g - second polynomial
    - var_elim - variable, of which f and g are considered as polynomials
Output:
    - M::MatrixElem, The Bezout matrix
"""
function Bezout_matrix(f::P, g::P, var_elim::P) where P <: MPolyElem
    parent_ring = parent(f)
    deg_f = degree(f, var_elim)
    deg_g = degree(g, var_elim)
    n = max(deg_f, deg_g)
    coeffs_f = [coeff(f, [var_elim], [i]) for i in 0:n]
    coeffs_g = [coeff(g, [var_elim], [i]) for i in 0:n]
    GL = AbstractAlgebra.MatrixSpace(parent_ring, n, n)
    M = zero(GL)
    for i in 1:n
        for j in 1:n
            M[i, j] = sum(coeffs_f[j + k + 1] * coeffs_g[i - k] - coeffs_g[j + k + 1] * coeffs_f[i - k] for k in 0:min(i - 1, n - j))
        end
    end
    return M
end

# ------------------------------------------------------------------------------

"""
    Sylvester_matrix(f, g, var_elim)

Compute the Sylvester matrix of two polynomials f, g with respect to var_elim
Inputs:
    - f - first polynomial
    - g - second polynomial
    - var_elim - variable, of which f and g are considered as polynomials
Output:
    - M::MatrixElem, The Sylvester matrix
"""
function Sylvester_matrix(f::P, g::P, var_elim::P) where P <: MPolyElem
    parent_ring = parent(f)
    deg_f = degree(f, var_elim)
    deg_g = degree(g, var_elim)
    n = deg_f + deg_g
    GL = AbstractAlgebra.MatrixSpace(parent_ring, n, n)
    M = zero(GL)
    for i in 1:deg_f
        for j in 0:deg_g
            M[i, j + i] = coeff(g, [var_elim], [j])
        end
    end
    for i in 1:deg_g
        for j in 0:deg_f
            M[i + deg_f, j + i] = coeff(f, [var_elim], [j])
        end
    end

    return M
end

# ------------------------------------------------------------------------------

"""
    simplify_matrix(M)

Eliminate GCD of entries of every row and column
Input:
    - M::MatrixElem, matrix to be simplified
Output:
    - M::MatrixElem, Simplified matrix
    - extra_factors::Vector{AbstractAlgebra.MPolyElem}, array of GCDs eliminated from M.
"""
function simplify_matrix(M::MatElem{P}) where P <: MPolyElem
    """
    An auxiliary function taking a list of coordinates of cells
    and dividing them by their gcd.
    Returns the gcd
    """
    function _simplify_range(coords::Array{Tuple{Int,Int},1})
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

    extra_factors = Array{P,1}()
    rows_cols = Array{Array{Tuple{Int,Int},1},1}()
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

# ------------------------------------------------------------------------------

# Definition: for an irreducible polynomial P in n variables, we will call 
# an iterator Julia object *generic point generator* if it generates 
# an infinite random sequence of points in C^n such that
#    1. P vanishes at each of the points
#    2. for every other irreducible polynomial, the probability of
#    vanishing at each of these points is zero
#
# Each returned point is dictionary from the variables to the values

# ------------------------------------------------------------------------------

abstract type PointGenerator{P} end

mutable struct ODEPointGenerator{P} <: PointGenerator{P}
    ode::ODE{P}
    big_ring::MPolyRing
    precision::Int
    cached_points::Array{Dict{P,<: FieldElem},1}
    number_type::Type

    function ODEPointGenerator{P}(ode::ODE{P}, big_ring::MPolyRing) where P <: MPolyElem
        prec = length(ode.x_vars) + 1
        number_type = typeof(one(base_ring(big_ring)))
        return new(ode, big_ring, prec, Array{Dict{P,number_type}}[], number_type)
    end
end

# ------------------------------------------------------------------------------

function Base.iterate(gpg::ODEPointGenerator{P}, i::Int=1) where P <: MPolyElem{<: FieldElem}
    if i > length(gpg.cached_points)
        @debug "Generating new point on the variety"
        sample_max = i * 50
        result = undef
        while true
            @debug "Preparing initial condition"
            flush(stdout)
            base_field = base_ring(gpg.big_ring)
            param_values = Dict{P,Int}(p => rand(1:sample_max) for p in gpg.ode.parameters)
            initial_conditions = Dict{P,Int}(x => rand(1:sample_max) for x in gpg.ode.x_vars)
            input_values = Dict{P,Array{Int,1}}(u => [rand(1:sample_max) for _ in 1:gpg.precision] for u in gpg.ode.u_vars)
            @debug "Computing a power series solution"
            flush(stdout)
            ps_solution = undef
            try
                ps_solution = power_series_solution(gpg.ode, param_values, initial_conditions, input_values, gpg.precision)
            catch e
                @debug "$e"
                flush(stdout)
                continue
            end

            @debug "Constructing the point"
            flush(stdout)
            result = Dict{P,gpg.number_type}(switch_ring(p, gpg.big_ring) => base_field(c) for (p, c) in param_values)
            for u in gpg.ode.u_vars
                for i in 0:(gpg.precision - 1)
                    result[str_to_var(var_to_str(u) * "_$i", gpg.big_ring)] = coeff(ps_solution[u], i) * factorial(i)
                end
            end
            for y in gpg.ode.y_vars
                for j in 0:(gpg.precision - 1)
                    result[str_to_var(var_to_str(y) * "_$j", gpg.big_ring)] = coeff(ps_solution[y], j) * factorial(j)
                end
            end
            for x in gpg.ode.x_vars
                result[switch_ring(x, gpg.big_ring)] = coeff(ps_solution[x], 0)
                result[str_to_var(var_to_str(x) * "_dot", gpg.big_ring)] = coeff(ps_solution[x], 1)
            end
            break
        end
        push!(gpg.cached_points, result)
    end
    return (gpg.cached_points[i], i + 1)
end


# ------------------------------------------------------------------------------

"""
    choose(polys, generic_point_generator)

Input:
    - polys, an array of distinct irreducible polynomials in the same ring
    - generic_point_generator, a generic point generator as described above for one of polys
Output:
    - the polynomial that vanishes at the generic_point_generator
"""
function choose(polys::Array{P,1}, generic_point_generator) where P <: MPolyElem{<: FieldElem}
    vars = gens(parent(polys[1]))
    for p in generic_point_generator
        if length(polys) <= 1
            break
        end
        # get accounts for the fact that the big ring may contain some auxiliary variables, e.g. rand_proj_var
        point = [get(p, v, zero(base_ring(parent(polys[1])))) for v in vars]
        polys = filter(e -> (evaluate(e, point) == 0), polys)
        flush(stdout)
    end
    return polys[1]
end

# ------------------------------------------------------------------------------

"""
    eliminate_var(f, g, var_elim, generic_point_generator)

Eliminate variable from a pair of polynomials
Input:
    - f and g, polynomials
    - var_elim, variable to be eliminated
    - generic_point_generator, a generic point generator object for the factor
      of the resultant of f and g of interest
Output:
    polynomial, the desired factor of the resultant of f and g
"""
function eliminate_var(f::P, g::P, var_elim::P, generic_point_generator) where P <: MPolyElem{<: FieldElem}
    # Linear comb
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
            g = g - q * f * var_elim^(degree(g, var_elim) - degree(f, var_elim))
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

    # Case (f(v^d), g(v^d)):
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
        f = sum(coeff(f, [var_elim], [gcd_deg * i]) * (var_elim^i) for i in 0:(degree(f, var_elim) ÷ gcd_deg))
        g = sum(coeff(g, [var_elim], [gcd_deg * i]) * (var_elim^i) for i in 0:(degree(g, var_elim) ÷ gcd_deg))
    end

    resultant = undef
    matrix_factors = []
    if degree(f, var_elim) == 0
        resultant = f
    else
        if degree(f, var_elim) > 1
            @debug "Calculating Bezout Matrix $(Dates.now())"
            flush(stdout)
            M = Bezout_matrix(f, g, var_elim)
        else
            @debug "Calculating Sylvester matrix"
            flush(stdout)
            M = Sylvester_matrix(f, g, var_elim)
        end
        @debug "Simplifying the matrix $(Dates.now())"
        flush(stdout)
        M_simp, matrix_factors = simplify_matrix(M)
        @debug "Removed factors $(map(length, matrix_factors))"
        M_size = zero(Nemo.MatrixSpace(Nemo.ZZ, ncols(M_simp), ncols(M_simp)))
        for i in 1:ncols(M_simp)
            for j in 1:ncols(M_simp)
                M_size[i,j] = length(M_simp[i,j])
            end
        end
        @debug "\t Matrix size: \n $M_size"
        @debug "\t Computing determinant $(Dates.now())"
        flush(stdout)
        resultant = det_minor_expansion(M_simp)
    end
    # Step 4: Eliminate extra factors
    factors = fast_factor(resultant)
    for mfac in matrix_factors
        for fac in fast_factor(mfac)
            if !(fac in factors)
                push!(factors, fac)
            end
        end
    end
    res = choose(factors, generic_point_generator)
    for f in factors
        if f != res
            @debug "\t \t Size of extra factor: $(length(f)); $(Dates.now())"
            @debug "\t \t It is $f"
            flush(stdout)
        end
    end
    return res
end
