"""
    sequence_solution(dds, param_values, initial_conditions, input_values, num_terms)

Input:
- `dds` - a discrete dynamical system to solve (represented as an ODE struct)
- `param_values` - parameter values, must be a dictionary mapping parameter to a value
- `initial_conditions` - initial conditions of `ode`, must be a dictionary mapping state variable to a value
- `input_values` - input sequences in the form input => list of terms; length of the lists must be at least
                   teh required number of terms in the result
- `num_terms` - number of terms to compute

Output:
- computes a sequence solution with teh required number of terms prec presented as a dictionary state_variable => corresponding sequence
"""
function sequence_solution(
    dds::ODE{P},
    param_values::Dict{P, T},
    initial_conditions::Dict{P, T},
    input_values::Dict{P, Array{T, 1}},
    num_terms::Int,
) where {T <: FieldElem, P <: MPolyElem{T}}
    result = Dict(x => [initial_conditions[x]] for x in dds.x_vars)
    for i in 2:num_terms
        eval_dict = merge(
            param_values, 
            Dict(k => v[end] for (k, v) in result),
            Dict(u => val[i - 1] for (u, val) in input_values)
        )
        for x in dds.x_vars
            push!(result[x], eval_at_dict(dds.x_equations[x], eval_dict))
        end
    end
    return result
end

#------------------------------------------------------------------------------

function sequence_solution(
    dds::ODE{P},
    param_values::Dict{P, Int},
    initial_conditions::Dict{P, Int},
    input_values::Dict{P, Array{Int, 1}},
    num_terms::Int,
) where {P <: MPolyElem{<:FieldElem}}
    bring = base_ring(dds.poly_ring)
    return sequence_solution(
        dds,
        Dict(p => bring(v) for (p, v) in param_values),
        Dict(x => bring(v) for (x, v) in initial_conditions),
        Dict(u => map(v -> bring(v), vv) for (u, vv) in input_values),
        num_terms
    )
end

#-----------------------------------------------------------------------------

"""
    differentiate_sequence_solution(dds, params, ic, inputs, num_terms)

Input:
- the same as for `sequence_solutions`

Output:
- a tuple consisting of the sequence solution and a dictionary of the form `(u, v) => sequence`, where `u` is a state variable
  `v` is a state or parameter, and the sequence is the partial derivative of
  the function `u` w.r.t. `v` evaluated at the solution
"""
function differentiate_sequence_solution(
    dds::ODE{P},
    params::Dict{P, T},
    ic::Dict{P, T},
    inputs::Dict{P, Array{T, 1}},
    num_terms::Int,
) where {T <: Generic.FieldElem, P <: MPolyElem{T}}
    @debug "Computing the power series solution of the system"
    seq_sol = sequence_solution(dds, params, ic, inputs, num_terms)
    generalized_params = vcat(dds.x_vars, dds.parameters)
    bring = base_ring(dds.poly_ring)

    @debug "Solving the variational system at the solution"
    part_diffs = Dict(
        (x, p) => derivative(dds.x_equations[x], p) 
        for x in dds.x_vars for p in generalized_params
    )
    result = Dict(
        (x, p) => [x == p ? one(bring) : zero(bring)]
        for x in dds.x_vars for p in generalized_params
    )
    for i in 2:num_terms
        eval_dict = merge(
            params, 
            Dict(k => v[i - 1] for (k, v) in seq_sol),
            Dict(u => val[i - 1] for (u, val) in inputs)
        )
        for p in generalized_params
            local_eval = Dict(x => result[(x, p)][end] for x in dds.x_vars)
            for x in dds.x_vars
                res = sum([eval_at_dict(part_diffs[(x, x2)], eval_dict) * local_eval[x2] for x2 in dds.x_vars])
                if p in dds.parameters
                    res += eval_at_dict(part_diffs[(x, p)], eval_dict)
                end
                push!(result[(x, p)], res)
            end
        end
    end

    return (seq_sol, result)
end

# ------------------------------------------------------------------------------

"""
    differentiate_sequence_output(dds, params, ic, inputs, num_terms)

Similar to `differentiate_sequence_solution` but computes partial derivatives of prescribed outputs
returns a dictionary of the form `y_function => Dict(var => dy/dvar)` where `dy/dvar` is the derivative
of `y_function` with respect to `var`.
"""
function differentiate_sequence_output(
    dds::ODE{P},
    params::Dict{P, T},
    ic::Dict{P, T},
    inputs::Dict{P, Array{T, 1}},
    num_terms::Int
) where {T <: Generic.FieldElem, P <: MPolyElem{T}}
    @debug "Computing partial derivatives of the solution"
    seq_sol, sol_diff = differentiate_sequence_solution(dds, params, ic, inputs, num_terms)

    @debug "Evaluating the partial derivatives of the outputs"
    generalized_params = vcat(dds.x_vars, dds.parameters)
    part_diffs = Dict(
        (y, p) => derivative(dds.y_equations[y], p) 
        for y in dds.y_vars for p in generalized_params
    )
 
    result = Dict((y, p) => [] for y in dds.y_vars for p in generalized_params)
    for i in 1:num_terms
        eval_dict = merge(
            params, 
            Dict(k => v[i] for (k, v) in seq_sol),
            Dict(u => val[i] for (u, val) in inputs)
        )
 
        for p in vcat(dds.x_vars, dds.parameters)
            local_eval = Dict(x => sol_diff[(x, p)][i] for x in dds.x_vars)
            for (y, y_eq) in dds.y_equations
                res = sum([eval_at_dict(part_diffs[(y, x)], eval_dict) * local_eval[x] for x in dds.x_vars])
                if p in dds.parameters
                    res += eval_at_dict(part_diffs[(y, p)], eval_dict)
                end
                push!(result[(y, p)], res)
            end
        end
    end

    return result
end

