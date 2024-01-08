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
            Dict(u => val[i - 1] for (u, val) in input_values),
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
        num_terms,
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
        (x, p) => derivative(dds.x_equations[x], p) for x in dds.x_vars for
        p in generalized_params
    )
    result = Dict(
        (x, p) => [x == p ? one(bring) : zero(bring)] for x in dds.x_vars for
        p in generalized_params
    )
    for i in 2:num_terms
        eval_dict = merge(
            params,
            Dict(k => v[i - 1] for (k, v) in seq_sol),
            Dict(u => val[i - 1] for (u, val) in inputs),
        )
        for p in generalized_params
            local_eval = Dict(x => result[(x, p)][end] for x in dds.x_vars)
            for x in dds.x_vars
                res = sum([
                    eval_at_dict(part_diffs[(x, x2)], eval_dict) * local_eval[x2] for
                    x2 in dds.x_vars
                ])
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
    num_terms::Int,
) where {T <: Generic.FieldElem, P <: MPolyElem{T}}
    @debug "Computing partial derivatives of the solution"
    seq_sol, sol_diff = differentiate_sequence_solution(dds, params, ic, inputs, num_terms)

    @debug "Evaluating the partial derivatives of the outputs"
    generalized_params = vcat(dds.x_vars, dds.parameters)
    part_diffs = Dict(
        (y, p) => derivative(dds.y_equations[y], p) for y in dds.y_vars for
        p in generalized_params
    )

    result = Dict((y, p) => [] for y in dds.y_vars for p in generalized_params)
    for i in 1:num_terms
        eval_dict = merge(
            params,
            Dict(k => v[i] for (k, v) in seq_sol),
            Dict(u => val[i] for (u, val) in inputs),
        )

        for p in vcat(dds.x_vars, dds.parameters)
            local_eval = Dict(x => sol_diff[(x, p)][i] for x in dds.x_vars)
            for (y, y_eq) in dds.y_equations
                res = sum([
                    eval_at_dict(part_diffs[(y, x)], eval_dict) * local_eval[x] for
                    x in dds.x_vars
                ])
                if p in dds.parameters
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
    _assess_local_identifiability_discrete_aux(dds::ODE{P}, funcs_to_check::Array{<: Any, 1}, known_ic, prob_threshold::Float64=0.99) where P <: MPolyElem{Nemo.fmpq}

Checks the local identifiability/observability of the functions in `funcs_to_check` treating `dds` as a discrete-time system with **shift**
instead of derivative in the right-hand side.
The result is correct with probability at least `prob_threshold`.
`known_ic` can take one of the following
 * `:none` - no initial conditions are assumed to be known
 * `:all` - all initial conditions are assumed to be known
 * a list of rational functions in states and parameters assumed to be known at t = 0
"""
function _assess_local_identifiability_discrete_aux(
    dds::ODE{P},
    funcs_to_check::Array{<:Any, 1},
    known_ic = :none,
    prob_threshold::Float64 = 0.99,
) where {P <: MPolyElem{Nemo.fmpq}}
    bring = base_ring(dds.poly_ring)

    @debug "Extending the model"
    dds_ext =
        add_outputs(dds, Dict("loc_aux_$i" => f for (i, f) in enumerate(funcs_to_check)))

    if known_ic == :none
        known_ic = []
    end
    if known_ic == :all
        known_ic = dds_ext.x_vars
    end

    @debug "Computing the observability matrix"
    prec = length(dds.x_vars) + length(dds.parameters)
    @debug "The truncation order is $prec"

    # Computing the bound from the Schwartz-Zippel-DeMilo-Lipton lemma
    deg_x = _degree_with_common_denom(values(dds.x_equations))
    deg_y = _degree_with_common_denom(values(dds.y_equations))
    deg_known = reduce(+, map(total_degree, known_ic), init = 0)
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
    params_vals = Dict(p => bring(rand(1:D)) for p in dds_ext.parameters)
    ic = Dict(x => bring(rand(1:D)) for x in dds_ext.x_vars)
    # TODO: parametric type instead of fmpq
    inputs = Dict{P, Array{fmpq, 1}}(
        u => [bring(rand(1:D)) for i in 1:prec] for u in dds_ext.u_vars
    )

    @debug "Computing the output derivatives"
    output_derivatives =
        differentiate_sequence_output(dds_ext, params_vals, ic, inputs, prec)

    @debug "Building the matrices"
    Jac = zero(
        Nemo.MatrixSpace(
            bring,
            length(dds.x_vars) + length(dds.parameters),
            1 + prec * length(dds.y_vars) + length(known_ic),
        ),
    )
    xs_params = vcat(dds_ext.x_vars, dds_ext.parameters)
    for (i, y) in enumerate(dds.y_vars)
        y = switch_ring(y, dds_ext.poly_ring)
        for j in 1:prec
            for (k, p) in enumerate(dds_ext.parameters)
                Jac[k, 1 + (i - 1) * prec + j] = output_derivatives[(y, p)][j]
            end
            for (k, x) in enumerate(dds_ext.x_vars)
                Jac[end - k + 1, 1 + (i - 1) * prec + j] = output_derivatives[(y, x)][j]
            end
        end
    end
    eval_point = merge(params_vals, ic)
    for (i, v) in enumerate(known_ic)
        for (k, p) in enumerate(dds_ext.parameters)
            Jac[k, end - i + 1] = eval_at_dict(derivative(v, p), eval_point)
        end
        for (k, x) in enumerate(dds_ext.x_vars)
            Jac[end - k + 1, end - i + 1] = eval_at_dict(derivative(v, x), eval_point)
        end
    end

    @debug "Computing the result"
    base_rank = LinearAlgebra.rank(Jac)
    result = OrderedDict{Any, Bool}()
    for i in 1:length(funcs_to_check)
        for (k, p) in enumerate(dds_ext.parameters)
            Jac[k, 1] =
                output_derivatives[(str_to_var("loc_aux_$i", dds_ext.poly_ring), p)][1]
        end
        for (k, x) in enumerate(dds_ext.x_vars)
            Jac[end - k + 1, 1] =
                output_derivatives[(str_to_var("loc_aux_$i", dds_ext.poly_ring), x)][1]
        end
        result[funcs_to_check[i]] = LinearAlgebra.rank(Jac) == base_rank
    end

    return result
end

# ------------------------------------------------------------------------------

"""
    function assess_local_identifiability(
        dds::ModelingToolkit.DiscreteSystem; 
        measured_quantities=Array{ModelingToolkit.Equation}[], 
        funcs_to_check=Array{}[], 
        known_ic=Array{}[],
        prob_threshold::Float64=0.99)

Input:
- `dds` - the DiscreteSystem object from ModelingToolkit (with **difference** operator in the right-hand side)
- `measured_quantities` - the measurable outputs of the model
- `funcs_to_check` - functions of parameters for which to check identifiability (all parameters and states if not specified)
- `known_ic` - functions (of states and parameter) whose initial conditions are assumed to be known
- `prob_threshold` - probability of correctness

Output:
- the result is an (ordered) dictionary from each function to to boolean;

The result is correct with probability at least `prob_threshold`.
"""
function assess_local_identifiability(
    dds::ModelingToolkit.DiscreteSystem;
    measured_quantities = Array{ModelingToolkit.Equation}[],
    funcs_to_check = Array{}[],
    known_ic = Array{}[],
    prob_threshold::Float64 = 0.99,
    loglevel = Logging.Info,
)
    restart_logging(loglevel = loglevel)
    with_logger(_si_logger[]) do
        return _assess_local_identifiability(
            dds,
            measured_quantities = measured_quantities,
            funcs_to_check = funcs_to_check,
            known_ic = known_ic,
            prob_threshold = prob_threshold,
        )
    end
end

function _assess_local_identifiability(
    dds::ModelingToolkit.DiscreteSystem;
    measured_quantities = Array{ModelingToolkit.Equation}[],
    funcs_to_check = Array{}[],
    known_ic = Array{}[],
    prob_threshold::Float64 = 0.99,
)
    if length(measured_quantities) == 0
        if any(ModelingToolkit.isoutput(eq.lhs) for eq in ModelingToolkit.equations(dds))
            @info "Measured quantities are not provided, trying to find the outputs in input dynamical system."
            measured_quantities = filter(
                eq -> (ModelingToolkit.isoutput(eq.lhs)),
                ModelingToolkit.equations(dds),
            )
        else
            throw(
                error(
                    "Measured quantities (output functions) were not provided and no outputs were found.",
                ),
            )
        end
    end

    # Converting the finite difference operator in the right-hand side to
    # the corresponding shift operator
    eqs = filter(eq -> !(ModelingToolkit.isoutput(eq.lhs)), ModelingToolkit.equations(dds))
    deltas = [Symbolics.operation(e.lhs).dt for e in eqs]
    @assert length(Set(deltas)) == 1
    eqs_shift = [e.lhs ~ e.rhs + first(Symbolics.arguments(e.lhs)) for e in eqs]
    dds_shift = DiscreteSystem(eqs_shift, name = gensym())
    @debug "System transformed from difference to shift: $dds_shift"

    dds_aux, conversion = mtk_to_si(dds_shift, measured_quantities)
    if length(funcs_to_check) == 0
        params = parameters(dds)
        params_from_measured_quantities = union(
            [filter(s -> !istree(s), get_variables(y)) for y in measured_quantities]...,
        )
        funcs_to_check = vcat(
            [x for x in states(dds) if conversion[x] in dds_aux.x_vars],
            union(params, params_from_measured_quantities),
        )
    end
    funcs_to_check_ = [eval_at_nemo(x, conversion) for x in funcs_to_check]
    known_ic_ = [eval_at_nemo(x, conversion) for x in known_ic]

    result = _assess_local_identifiability_discrete_aux(
        dds_aux,
        funcs_to_check_,
        known_ic_,
        prob_threshold,
    )
    nemo2mtk = Dict(funcs_to_check_ .=> funcs_to_check)
    out_dict = OrderedDict(nemo2mtk[param] => result[param] for param in funcs_to_check_)
    if length(known_ic) > 0
        @warn "Since known initial conditions were provided, identifiability of states (e.g., `x(t)`) is at t = 0 only !"
        out_dict = OrderedDict(substitute(k, Dict(t => 0)) => v for (k, v) in out_dict)
    end
    return out_dict
end
# ------------------------------------------------------------------------------
