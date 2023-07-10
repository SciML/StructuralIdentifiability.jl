function extract_coefficients_ratfunc(f::AbstractAlgebra.Generic.Frac{<: P}, vars::Vector{<: P}) where {P <: MPolyElem{<:FieldElem}}
    num, denom = unpack_fraction(f)
    total_coeffs = Vector{P}()
    for p in (num, denom)
        coeffs = map(c -> parent_ring_change(c, parent(p)), values(extract_coefficients(p, vars)))
        append!(total_coeffs, coeffs)
    end
    divisor = total_coeffs[argmin(map(c -> (total_degree(c), length(c)), total_coeffs))]
    return [x // divisor for x in total_coeffs]
end
