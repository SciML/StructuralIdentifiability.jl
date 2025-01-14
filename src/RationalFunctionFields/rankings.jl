# Rank field generators by their simplicity

"""
    rational_function_rank(func)

The rank of the given rational function, a nonnegative integer. 
Smaller rank is better.
"""
function rational_function_rank(func)
    num, den = numerator(func), denominator(func)
    r1 = length(num) * total_degree(num)
    # penalize non-constant denominators
    r2 = (3 * length(den) * total_degree(den))^3
    # penalize coefficients that are not 1
    r3 = sum(map(c -> !isone(c), collect(coefficients(num))))
    r4 = sum(map(c -> !isone(c), collect(coefficients(den))))
    return r1 + r2 + r3 + r4
end

"""
    generating_set_rank(funcs)

The rank of the given set of generators, a nonnegative integer. 
Smaller rank is better.
"""
function generating_set_rank(funcs)
    n = length(funcs)
    func_ranks = Vector{Int}(undef, n)
    for i in 1:n
        func_ranks[i] = rational_function_rank(funcs[i])
    end
    rank = BigInt(0)
    for i in 1:n
        rank = rank + (n - i + 1) * func_ranks[i]
    end
    return rank
end

"""
    rational_function_cmp(f, g; by=:naive, priority_variables=[])

Returns `true` if `f < g`.

*This is a strict total order, which ensures the uniqueness of the sorting
permutation.*

Provides keyword argument `by`, a sorting order. Possible options are:
- `:rank`: Compare by rank.
- `:naive`: Compare features one by one. Features in the order of importance:
    - Sum of total degrees
    - Total number of terms
    - Total degree of the denominator
    - Leading terms of denominator and numerator
"""
function rational_function_cmp(f, g; by = :naive)
    if by === :naive
        flag = compare_rational_func_by(f, g, total_degree, :additive)
        flag == 1 && return false
        flag == -1 && return true
        flag = compare_rational_func_by(f, g, length, :additive)
        flag == 1 && return false
        flag == -1 && return true
        flag = compare_rational_func_by(f, g, total_degree, :denominator)
        flag == 1 && return false
        flag == -1 && return true
        # promotes constants in denominators
        flag = compare_rational_func_by(f, g, leading_monomial, :denominator)
        flag == 1 && return false
        flag == -1 && return true
        flag = compare_rational_func_by(f, g, collect ∘ monomials)
        flag == 1 && return false
        flag == -1 && return true
        return false
    else
        @assert by === :rank
        return rational_function_rank(f) < rational_function_rank(g)
    end
end
