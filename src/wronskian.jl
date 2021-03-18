import Base.push!

#----------------------------------------------------------------------------------------------------

"""
    monomial_compress(io_equation, ode)

Compresses an input-output equation for the rank computation
Input: 
    - io_equation - input-output equation
    - ode - the corresponding ode model
Output: pair (coeffs, terms) such that
    - sum of coeffs[i] * terms[i] = io_equation
    - coeffs involve only parameters, terms involve only inputs and outputs
    - length of the representation is the smallest possible
"""
function monomial_compress(io_equation, ode::ODE)
	params = ode.parameters
    other_vars = [v for v in gens(parent(io_equation)) if !(var_to_str(v) in map(var_to_str, params))]
	coeffdict = extract_coefficients(io_equation, other_vars)
	expvect = collect(keys(coeffdict))
	coeffs = collect(values(coeffdict))
	termlist = map(x->prod(other_vars.^x), expvect)

    echelon_form = Dict()
    for (c, p) in zip(coeffs, termlist)
        for basis_c in keys(echelon_form)
            coef = coeff(c, lm(basis_c)) // lc(basis_c)
            if coef != 0
                c = c - coef * basis_c
                echelon_form[basis_c] += coef * p
            end
        end
        if c != 0
            echelon_form[c] = p
        end
    end

    return (collect(keys(echelon_form)), collect(values(echelon_form)))
end

#----------------------------------------------------------------------------------------------------

# A data structure for storing exponent vectors as a trie in order
# to reduce the number of multiplications when computing many monomials
# in `massive_eval`
mutable struct ExpVectTrie
    depth::Int
    subtries::Dict{Int, ExpVectTrie}

    function ExpVectTrie(depth::Int)
        return new(depth, Dict{Int, ExpVectTrie}())
    end
end

function Base.push!(t::ExpVectTrie, vect::Array{Int, 1})
    if t.depth == 0
        if length(vect) != 0
            throw(DomainError("Inserting too long vector"))
        end
        return
    end

    f = last(vect)
    if !(f in keys(t.subtries))
        t.subtries[f] = ExpVectTrie(t.depth - 1)
    end
    push!(t.subtries[f], vect[1:end - 1])
end

"""
    get_max_below(t, vect)

Input:
    - t - a trie with exponent vectors
    - vect -yet another exponent vector
Output: a pair d, v, where v is a vector in the trie which is
    componenwise <= vect and the difference d is as small as possible
"""
function get_max_below(t::ExpVectTrie, vect::Array{Int, 1})
    if t.depth == 0
        return (0, [])
    end
    min_diff, best_below = nothing, nothing
    for k in keys(t.subtries)
        last_diff = last(vect) - k
        if last_diff >= 0
            cur_diff, cur_below = get_max_below(t.subtries[k], vect[1:end - 1])
            if !isnothing(cur_diff)
                if isnothing(min_diff) || min_diff > cur_diff + last_diff
                    min_diff = cur_diff + last_diff
                    push!(cur_below, k)
                    best_below = cur_below
                end
            end
        end
    end
    return (min_diff, best_below)
end

#----------------------------------------------------------------------------------------------------

"""
    massive_eval(polys, eval_dict)

Evaluates a list of polynomails at a point. Assumes that multiplications are relatively expensive
(like in truncated power series) so all the monomials are precomputed first and the values of monomials
of lower degree are cached and used to compute the values of the monomials of higher degree
Input:
    - polys - a list of polynomials
    - eval_dict - dictionary from variables to the values. Missing are treated as zeroes
Output: a list of values of the polynomials
"""
function massive_eval(polys, eval_dict)
    R = parent(first(values(eval_dict)))
    point = [get(eval_dict, v, zero(R)) for v in gens(parent(first(polys)))]
    n = length(point)

    monomials = Set()
    for p in polys
        for exp in exponent_vectors(p)
            push!(monomials, exp)
        end
    end

    cache = Dict()
    cache[ [0 for i in 1:n] ] = one(R)
    cached_monoms = ExpVectTrie(n)
    push!(cached_monoms, [0 for _ in 1:n])

    for i in 1:n
        var_exp = [(i != j) ? 0 : 1 for j in 1:n]
        cache[var_exp] = point[i]
        push!(cached_monoms, var_exp)
    end

    for exp in sort!(collect(monomials), by=sum)
        if !(exp in keys(cache))
            monom_val = one(R)
            computed = [0 for i in 1:n]
            while sum(exp) > 0
                _, below = get_max_below(cached_monoms, exp)
                monom_val = monom_val * cache[below]
                exp = exp .- below
                computed = computed .+ below
                cache[computed] = monom_val
                push!(cached_monoms, computed)
            end
        end
    end

    results = []
    for p in polys
        res = zero(R)
        for (exp, coef) in zip(exponent_vectors(p), coefficients(p))
            res += coef * cache[exp]
        end
        push!(results, res)
    end
    return results
end

#----------------------------------------------------------------------------------------------------

"""
    wronskian(io_equations, ode)

Computes the wronskians of io_equations
Input:
    - io_equations - a set of io-equations in the form of the Dict as returned by find_ioequations
    - ode - the ode object
Output: a list of wronskians evaluated at a point modulo prime
"""
function wronskian(io_equations::Dict{P, P}, ode::ODE{P}) where P <: MPolyElem
    @debug "Compressing monomials"
    termlists = [monomial_compress(ioeq, ode)[2] for ioeq in values(io_equations)]
    @debug "Matrix sizes $(map(length, termlists))"

    # estimating the required order of truncation
    ord = max(map(length, termlists)...) + length(ode.x_vars)

    # reducing everything modulo prime
    PRIME = 2^31 - 1
    F = Nemo.GF(PRIME)
    polyring_red, gens_red = PolynomialRing(F, map(var_to_str, gens(parent(termlists[1][1]))))
    termlists = [map(p -> parent_ring_change(p, polyring_red), tlist) for tlist in termlists]
    ode_red = reduce_ode_mod_p(ode, PRIME) 

    @debug "Computing power series solution up to order $ord"
	ps = power_series_solution(
        ode_red,
        Dict(p => rand(1:100) for p in ode_red.parameters),
        Dict(x => rand(1:100) for x in ode_red.x_vars),
        Dict(u => [rand(1:100) for i in 0:ord] for u in ode_red.u_vars),
        ord
    )
    @debug "Computing the derivatives of the solution"
    ps_ext = Dict{MPolyElem, Generic.AbsSeries}()
    for v in vcat(ode_red.y_vars, ode_red.u_vars)
        cur_ps = ps[v]
        for i in 0:length(ode_red.x_vars)
            ps_ext[str_to_var(var_to_str(v) * "_$i", polyring_red)] = cur_ps
            cur_ps = ps_diff(cur_ps)
        end
    end

    @debug "Constructing wronskians"
    result = []
    for (i, tlist) in enumerate(termlists)
        n = length(tlist)
        evaled = massive_eval(tlist, ps_ext)
        S = Nemo.MatrixSpace(F, n, n)
        W = zero(S)
        for (i, ps) in enumerate(evaled)
            for j in 1:n
                W[i, j] = coeff(ps, j - 1)
            end
        end
        push!(result, W)
    end

    return result
end
