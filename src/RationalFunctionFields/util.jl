
"""
    dennums_to_fractions(dennums)
    
Returns the field generators represented by fractions.

Input: an array of arrays of polynomials, as in 
`[[f1, f2, f3, ...], [g1, g2, g3, ...], ...]`

Output: an array of fractions
`[f2/f1, f3/f1, ..., g2/g1, g3/g1, ...]`
"""
function dennums_to_fractions(dennums::Vector{Vector{T}}) where {T}
    fractions = Vector{AbstractAlgebra.Generic.Frac{T}}()
    for dni in dennums
        den, nums = dni[1], dni[2:end]
        isempty(nums) && continue
        append!(fractions, map(c -> c // den, nums))
    end
    return fractions
end

# ------------------------------------------------------------------------------

"""
    fractions_to_dennums(fractions)
    
Returns the field generators represented by lists of denominators and
numerators.

Input: an array of fractions, as in
`[f2/f1, f3/f1, ..., g2/g1, g3/g1, ...]`

Output: an array of arrays of polynomials,
`[[f1, f2, f3, ...], [g1, g2, g3, ...], ...]`
"""
function fractions_to_dennums(fractions)
    return map(f -> [denominator(f), numerator(f)], fractions)
end
