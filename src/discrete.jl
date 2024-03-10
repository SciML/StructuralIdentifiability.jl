"""
The structure to represent a discrete dynamical system
with respect to *shift*. Internally just stores an ODE structure

Can be constructed with @DDSmodel macro
"""
struct DDS{P} # P is the type of polynomials in the rhs of the DDS system
    ode::ODE{P}

    function DDS{P}(
        x_vars::Array{P, 1},
        y_vars::Array{P, 1},
        x_eqs::Dict{P, <:ExtendedFraction{P}},
        y_eqs::Dict{P, <:ExtendedFraction{P}},
        u_vars::Array{P, 1},
    ) where {P <: MPolyRingElem{<:FieldElem}}
        new{P}(ODE{P}(x_vars, y_vars, x_eqs, y_eqs, u_vars))
    end

    function DDS{P}(ode::ODE{P}) where {P <: MPolyRingElem{<:FieldElem}}
        new{P}(ode)
    end
end

#------------------------------------------------------------------------------

# getters

function x_vars(dds::DDS)
    return dds.ode.x_vars
end

function y_vars(dds::DDS)
    return dds.ode.y_vars
end

function parameters(dds::DDS)
    return dds.ode.parameters
end

function inputs(dds::DDS)
    return dds.ode.u_vars
end

function x_equations(dds::DDS)
    return dds.ode.x_equations
end

function y_equations(dds::DDS)
    return dds.ode.y_equations
end

function Base.parent(dds::DDS)
    return parent(dds.ode)
end

#------------------------------------------------------------------------------
# Some functions to transform DDS's

function add_outputs(
    dds::DDS{P},
    extra_y::Dict{String, <:RingElem},
) where {P <: MPolyRingElem}
    return DDS{P}(add_outputs(dds.ode, extra_y))
end

#------------------------------------------------------------------------------

function Base.show(io::IO, dds::DDS)
    for x in x_vars(dds)
        if endswith(var_to_str(x), "(t)")
            print(io, chop(var_to_str(x), tail = 3) * "(t + 1) = ")
        else
            print(io, var_to_str(x) * "(t + 1) = ")
        end
        print(io, x_equations(dds)[x])
        print(io, "\n")
    end
    for y in y_vars(dds)
        print(io, var_to_str(y) * " = ")
        print(io, y_equations(dds)[y])
        print(io, "\n")
    end
end

#------------------------------------------------------------------------------

"""
    sequence_solution(dds, param_values, initial_conditions, input_values, num_terms)

Input:
- `dds` - a discrete dynamical system to solve
- `param_values` - parameter values, must be a dictionary mapping parameter to a value
- `initial_conditions` - initial conditions of `ode`, must be a dictionary mapping state variable to a value
- `input_values` - input sequences in the form input => list of terms; length of the lists must be at least
                   the required number of terms in the result
- `num_terms` - number of terms to compute

Output:
- computes a sequence solution with the required number of terms prec presented as a dictionary state_variable => corresponding sequence
"""
function sequence_solution(
    dds::DDS{P},
    param_values::Dict{P, T},
    initial_conditions::Dict{P, T},
    input_values::Dict{P, Array{T, 1}},
    num_terms::Int,
) where {T <: FieldElem, P <: MPolyRingElem{T}}
    result = Dict(x => [initial_conditions[x]] for x in x_vars(dds))
    for i in 2:num_terms
        eval_dict = merge(
            param_values,
            Dict(k => v[end] for (k, v) in result),
            Dict(u => val[i - 1] for (u, val) in input_values),
        )
        for x in x_vars(dds)
            push!(result[x], eval_at_dict(x_equations(dds)[x], eval_dict))
        end
    end
    return result
end

#------------------------------------------------------------------------------

function sequence_solution(
    dds::DDS{P},
    param_values::Dict{P, Int},
    initial_conditions::Dict{P, Int},
    input_values::Dict{P, Array{Int, 1}},
    num_terms::Int,
) where {P <: MPolyRingElem{<:FieldElem}}
    bring = base_ring(parent(dds))
    return sequence_solution(
        dds,
        Dict(p => bring(v) for (p, v) in param_values),
        Dict(x => bring(v) for (x, v) in initial_conditions),
        Dict(u => map(v -> bring(v), vv) for (u, vv) in input_values),
        num_terms,
    )
end

#-----------------------------------------------------------------------------

"""
    differentiate_sequence_solution(dds, params, ic, input_values, num_terms)

Input:
- the same as for `sequence_solutions`

Output:
- a tuple consisting of the sequence solution and a dictionary of the form `(u, v) => sequence`, where `u` is a state variable
  `v` is a state or parameter, and the sequence is the partial derivative of
  the function `u` w.r.t. `v` evaluated at the solution
"""
function differentiate_sequence_solution(
    dds::DDS{P},
    params::Dict{P, T},
    ic::Dict{P, T},
    input_values::Dict{P, Array{T, 1}},
    num_terms::Int,
) where {T <: Generic.FieldElem, P <: MPolyRingElem{T}}
    @debug "Computing the power series solution of the system"
    seq_sol = sequence_solution(dds, params, ic, input_values, num_terms)
    generalized_params = vcat(x_vars(dds), parameters(dds))
    bring = base_ring(parent(dds))

    @debug "Solving the variational system at the solution"
    part_diffs = Dict(
        (x, p) => derivative(x_equations(dds)[x], p) for x in x_vars(dds) for
        p in generalized_params
    )
    result = Dict(
        (x, p) => [x == p ? one(bring) : zero(bring)] for x in x_vars(dds) for
        p in generalized_params
    )
    for i in 2:num_terms
        eval_dict = merge(
            params,
            Dict(k => v[i - 1] for (k, v) in seq_sol),
            Dict(u => val[i - 1] for (u, val) in input_values),
        )
        for p in generalized_params
            local_eval = Dict(x => result[(x, p)][end] for x in x_vars(dds))
            for x in x_vars(dds)
                res = sum([
                    eval_at_dict(part_diffs[(x, x2)], eval_dict) * local_eval[x2] for
                    x2 in x_vars(dds)
                ])
                if p in parameters(dds)
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
    differentiate_sequence_output(dds, params, ic, input_values, num_terms)

Similar to `differentiate_sequence_solution` but computes partial derivatives of prescribed outputs
returns a dictionary of the form `y_function => Dict(var => dy/dvar)` where `dy/dvar` is the derivative
of `y_function` with respect to `var`.
"""
function differentiate_sequence_output(
    dds::DDS{P},
    params::Dict{P, T},
    ic::Dict{P, T},
    input_values::Dict{P, Array{T, 1}},
    num_terms::Int,
) where {T <: Generic.FieldElem, P <: MPolyRingElem{T}}
    @debug "Computing partial derivatives of the solution"
    seq_sol, sol_diff =
        differentiate_sequence_solution(dds, params, ic, input_values, num_terms)

    @debug "Evaluating the partial derivatives of the outputs"
    generalized_params = vcat(x_vars(dds), parameters(dds))
    part_diffs = Dict(
        (y, p) => derivative(y_equations(dds)[y], p) for y in y_vars(dds) for
        p in generalized_params
    )

    result = Dict((y, p) => [] for y in y_vars(dds) for p in generalized_params)
    for i in 1:num_terms
        eval_dict = merge(
            params,
            Dict(k => v[i] for (k, v) in seq_sol),
            Dict(u => val[i] for (u, val) in input_values),
        )

        for p in generalized_params
            local_eval = Dict(x => sol_diff[(x, p)][i] for x in x_vars(dds))
            for (y, y_eq) in y_equations(dds)
                res = sum([
                    eval_at_dict(part_diffs[(y, x)], eval_dict) * local_eval[x] for
                    x in x_vars(dds)
                ])
                if p in parameters(dds)
                    res += eval_at_dict(part_diffs[(y, p)], eval_dict)
                end
                push!(result[(y, p)], res)
            end
        end
    end

    return result
end

# ------------------------------------------------------------------------------

function _degree_with_common_denom(polys)
    common_denom = reduce(lcm, map(f -> unpack_fraction(f)[2], polys))
    return max(
        max(
            [
                total_degree(pair[1]) + total_degree(common_denom) - total_degree(pair[2])
                for pair in map(unpack_fraction, polys)
            ]...,
        ),
        total_degree(common_denom),
    )
end

"""
    _assess_local_identifiability_discrete_aux(dds::DDS{P}, funcs_to_check::Array{<: Any, 1}, known_ic, prob_threshold::Float64=0.99) where P <: MPolyRingElem{Nemo.QQFieldElem}

Checks the local identifiability/observability of the functions in `funcs_to_check`.
The result is correct with probability at least `prob_threshold`.
`known_ic` can take one of the following
 * `:none` - no initial conditions are assumed to be known
 * `:all` - all initial conditions are assumed to be known
 * a list of rational functions in states and parameters assumed to be known at t = 0
"""
function _assess_local_identifiability_discrete_aux(
    dds::DDS{P},
    funcs_to_check::Array{<:Any, 1},
    known_ic = :none,
    prob_threshold::Float64 = 0.99,
) where {P <: MPolyRingElem{Nemo.QQFieldElem}}
    bring = base_ring(parent(dds))

    @debug "Extending the model"
    dds_ext =
        add_outputs(dds, Dict("loc_aux_$i" => f for (i, f) in enumerate(funcs_to_check)))

    if known_ic == :none
        known_ic = []
    end
    if known_ic == :all
        known_ic = x_vars(dds_ext)
    end

    @debug "Computing the observability matrix"
    prec = length(x_vars(dds)) + length(parameters(dds))
    @debug "The truncation order is $prec"

    # Computing the bound from the Schwartz-Zippel-DeMilo-Lipton lemma
    deg_x = _degree_with_common_denom(values(x_equations(dds)))
    deg_y = _degree_with_common_denom(values(y_equations(dds)))
    deg_known = reduce(+, map(total_degree_frac, known_ic), init = 0)
    deg_to_check = max(map(total_degree_frac, funcs_to_check)...)
    Jac_degree = deg_to_check + deg_known
    if deg_x > 1
        Jac_degree += 2 * deg_y * ((deg_x^prec - 1) ÷ (deg_x - 1))
    else
        Jac_degree += 2 * deg_y * prec
    end
    D = Int(ceil(Jac_degree * length(funcs_to_check) / (1 - prob_threshold)))
    @debug "Sampling range $D"

    # Parameter values are the same across all the replicas
    params_vals = Dict(p => bring(rand(1:D)) for p in parameters(dds_ext))
    ic = Dict(x => bring(rand(1:D)) for x in x_vars(dds_ext))
    # TODO: parametric type instead of QQFieldElem
    input_values = Dict{P, Array{QQFieldElem, 1}}(
        u => [bring(rand(1:D)) for i in 1:prec] for
        u in StructuralIdentifiability.inputs(dds_ext)
    )

    @debug "Computing the output derivatives"
    output_derivatives =
        differentiate_sequence_output(dds_ext, params_vals, ic, input_values, prec)

    @debug "Building the matrices"
    Jac = zero(
        Nemo.matrix_space(
            bring,
            length(x_vars(dds)) + length(parameters(dds)),
            1 + prec * length(y_vars(dds)) + length(known_ic),
        ),
    )
    xs_params = vcat(x_vars(dds_ext), parameters(dds_ext))
    for (i, y) in enumerate(y_vars(dds))
        y = switch_ring(y, parent(dds_ext))
        for j in 1:prec
            for (k, p) in enumerate(parameters(dds_ext))
                Jac[k, 1 + (i - 1) * prec + j] = output_derivatives[(y, p)][j]
            end
            for (k, x) in enumerate(x_vars(dds_ext))
                Jac[end - k + 1, 1 + (i - 1) * prec + j] = output_derivatives[(y, x)][j]
            end
        end
    end
    eval_point = merge(params_vals, ic)
    for (i, v) in enumerate(known_ic)
        for (k, p) in enumerate(parameters(dds_ext))
            Jac[k, end - i + 1] = eval_at_dict(derivative(v, p), eval_point)
        end
        for (k, x) in enumerate(x_vars(dds_ext))
            Jac[end - k + 1, end - i + 1] = eval_at_dict(derivative(v, x), eval_point)
        end
    end

    @debug "Computing the result"
    base_rank = LinearAlgebra.rank(Jac)
    result = OrderedDict{Any, Bool}()
    for i in 1:length(funcs_to_check)
        for (k, p) in enumerate(parameters(dds_ext))
            Jac[k, 1] =
                output_derivatives[(str_to_var("loc_aux_$i", parent(dds_ext)), p)][1]
        end
        for (k, x) in enumerate(x_vars(dds_ext))
            Jac[end - k + 1, 1] =
                output_derivatives[(str_to_var("loc_aux_$i", parent(dds_ext)), x)][1]
        end
        result[funcs_to_check[i]] = LinearAlgebra.rank(Jac) == base_rank
    end

    return result
end

# ------------------------------------------------------------------------------

"""
    assess_local_identifiability(dds::DDS{P}; funcs_to_check::Array{<: Any, 1}, known_ic, prob_threshold::Float64=0.99, loglevel=Logging.Info) where P <: MPolyRingElem{Nemo.QQFieldElem}

Checks the local identifiability/observability of the functions in `funcs_to_check`. The result is correct with probability at least `prob_threshold`.
A list of quantities can be provided as `known_ic` for which the initial conditions can be assumed to be known and generic.
"""
function assess_local_identifiability(
    dds::DDS{P};
    funcs_to_check::Array{<:Any, 1} = Array{Any, 1}(),
    known_ic = :none,
    prob_threshold::Float64 = 0.99,
    loglevel = Logging.Info,
) where {P <: MPolyRingElem{Nemo.QQFieldElem}}
    restart_logging(loglevel = loglevel)
    reset_timings()
    with_logger(_si_logger[]) do
        if isempty(funcs_to_check)
            funcs_to_check = vcat(parameters(dds), x_vars(dds))
        end
        return _assess_local_identifiability_discrete_aux(
            dds,
            funcs_to_check,
            known_ic,
            prob_threshold,
        )
    end
end

# ------------------------------------------------------------------------------
