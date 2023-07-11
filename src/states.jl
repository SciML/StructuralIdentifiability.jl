"""
    extract_coefficients_ratfunc(f, vars)

Input:
- `f` - rational function
- `vars` - list of variables

The function considers `f` as `A / B`, where `A` and `B` are polynomials in `vars` with
coefficients in rational fucntion field in the remaining variables such that at least one of the
coefficients is equal to one.

Output:
- list containing coefficients of `A` and `B` in some order (and in the original ring!)
"""
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

#------------------------------------------------------------------------------
"""
    lie_derivative(f, ode)

Input:
- `f` - rational function in states, parameters  (not inputs) of `ode`
- `ode' - an ODE model

Output:
- Lie derivative of `f` with respect to `ode` 
"""
function lie_derivative(f::Generic.Frac{<: P}, ode::ODE{<: P}) where {P <: MPolyElem{<: FieldElem}}
    @assert all([(x in ode.parameters) || (x in ode.x_vars) for x in vars(f)])
    res = zero(parent(ode)) // 1
    for (x, eq) in ode.x_equations
        res += derivative(f, x) * eq
    end
    return res
end

function lie_derivative(f::P, ode::ODE{<: P}) where {P <: MPolyElem{<: FieldElem}}
    return lie_derivative(f // 1, ode)
end
