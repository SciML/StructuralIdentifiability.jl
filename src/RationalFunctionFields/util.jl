
"""
    dennums_to_fractions(dennums)
    
Returns the field generators represented by fractions.

Input: an array of arrays of polynomials, as in 
`[[f1, f2, f3, ...], [g1, g2, g3, ...], ...]`

Output: an array of fractions
`[f2/f1, f3/f1, ..., g2/g1, g3/g1, ...]`
"""
function dennums_to_fractions(dennums::Vector{Vector{T}}) where {T}
    fractions = Vector{AbstractAlgebra.Generic.FracFieldElem{T}}()
    for dni in dennums
        pivot_ind = findmin(p -> (total_degree(p), length(p)), filter(!iszero, dni))[2]
        pivot = dni[pivot_ind]
        append!(fractions, [c // pivot for c in dni if c != pivot])
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

# ------------------------------------------------------------------------------

"""
    merge_results(outer::Vector{Bool}, inner::Vector{Bool})

Returns a list `res` of Bools of the length as outer such that `res[i]` is true
iff `outer[i]` is true and `inner[j]` is true where `j` is the index of `outer[i]`
among the true values in `outer`
"""

function merge_results(outer::Vector{Bool}, inner::Vector{Bool})
    @assert length(inner) == count(outer)
    result = copy(outer)
    inner_index = 1

    for i in 1:length(result)
        if result[i]
            if !inner[inner_index]
                result[i] = false
            end
            inner_index += 1
        end
    end
    return result
end
