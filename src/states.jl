"""
    extract_coefficients_ratfunc(f, vars)

Input:
- `f` - rational function
- `vars` - list of variables

The function considers `f` as `A / B`, where `A` and `B` are polynomials in `vars` with
coefficients in rational function field in the remaining variables such that at least one of the
coefficients is equal to one.

Output:
- list containing coefficients of `A` and `B` in some order (and in the original ring!)
"""
function extract_coefficients_ratfunc(
    f::AbstractAlgebra.Generic.FracFieldElem{<:P},
    vars::Vector{<:P},
) where {P <: MPolyRingElem{<:FieldElem}}
    num, denom = unpack_fraction(f)
    total_coeffs = Vector{P}()
    for p in (num, denom)
        coeffs = map(
            c -> parent_ring_change(c, parent(p)),
            values(extract_coefficients(p, vars)),
        )
        append!(total_coeffs, coeffs)
    end
    divisor = total_coeffs[argmin(map(c -> (total_degree(c), length(c)), total_coeffs))]
    return [x // divisor for x in total_coeffs]
end

function extract_coefficients_ratfunc(
    f::P,
    vars::Vector{<:P},
) where {P <: MPolyRingElem{<:FieldElem}}
    return extract_coefficients_ratfunc(f // 1, vars)
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
function lie_derivative(
    f::Generic.FracFieldElem{<:P},
    ode::ODE{<:P},
) where {P <: MPolyRingElem{<:FieldElem}}
    @assert all([(x in ode.parameters) || (x in ode.x_vars) for x in vars(f)])
    res = zero(parent(ode)) // 1
    for (x, eq) in ode.x_equations
        res += derivative(f, x) * eq
    end
    return res
end

function lie_derivative(f::P, ode::ODE{<:P}) where {P <: MPolyRingElem{<:FieldElem}}
    return lie_derivative(f // 1, ode)
end

#------------------------------------------------------------------------------
"""
    states_generators(ode, io_equations)

Input:
- `ode' - an ODE model
- `io_equations` - input-output equations of the model

Output:
- a list of rational functions in parameters and states generating the field of
all identifiabile functions over the field generated by the inputs and the
identifiable functions of parameters only
"""
@timeit _to function states_generators(
    ode::ODE{P},
    io_equations::Dict{P, P},
) where {P <: MPolyRingElem{<:FieldElem}}
    y_to_ord = Dict{P, Int}()
    ynames = [var_to_str(y) for y in ode.y_vars]
    for (leader, ioeq) in io_equations
        decomposition = decompose_derivative(var_to_str(leader), ynames)
        # decomposition will be `nothing` for extra random projection
        if !isnothing(decomposition)
            (y_str, ord) = decomposition
            y_to_ord[str_to_var(y_str, parent(ode))] = ord
        end
    end

    result = Array{Generic.FracFieldElem{P}, 1}()
    for (y, ord) in y_to_ord
        curr = extract_coefficients_ratfunc(ode.y_equations[y], ode.u_vars)
        for _ in 0:ord
            filter!(!is_rational_func_const, curr)
            append!(result, curr)
            curr = reduce(
                vcat,
                [
                    extract_coefficients_ratfunc(lie_derivative(f, ode), ode.u_vars) for
                    f in curr
                ],
                init = empty(curr),
            )
        end
    end

    return result
end
