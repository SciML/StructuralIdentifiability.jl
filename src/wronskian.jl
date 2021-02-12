using Dates
using IterTools
using Logging
using Oscar
using Base.Iterators

include("util.jl")
include("elimination.jl")
include("ODE.jl")
include("io_equation.jl")

#----------------------------------------------------------------------------------------------------

function monomial_compress(io_equation, ode::ODE)
    """
    Compresses an input-output equation for the rank computation
    Input: 
        - io_equation - input-output equation
        - ode - the corresponding ode model
    Output: pair (coeffs, terms) such that
        - sum of coeffs[i] * terms[i] = io_equation
        - coeffs involve only parameters, terms involve only inputs and outputs
        - length of the representation is the smallest possible
    """
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

# TODO generator
function massive_eval(polys, eval_dict)
    R = parent(first(values(eval_dict)))
    point = [get(eval_dict, v, zero(R)) for v in gens(parent(first(polys)))]

    cache = Dict()
    cache[ [0 for i in 1:length(point)] ] = one(R)

    results = []
    for p in polys
        res = zero(R)
        for (exp, coef) in zip(exponent_vectors(p), coefficients(p))
            if !(exp in keys(cache))
                prod = one(R)
                for (i, pow) in enumerate(exp)
                    if pow > 0
                        prod *= point[i]^pow
                    end
                end
                cache[exp] = prod
            end
            res += coef * cache[exp]
        end
        push!(results, res)
    end
    return results
end

#----------------------------------------------------------------------------------------------------

function wronskian(io_equations::Dict{P, P}, ode::ODE{P}) where P <: MPolyElem
    """
    Computes the wronskians of io_equations
    Input:
        - io_equations - a set of io-equations in the form of the Dict as returned by find_ioequations
        - ode - the ode object
    Output: a list of wronskians evaluated at a point modulo prime
    """
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
