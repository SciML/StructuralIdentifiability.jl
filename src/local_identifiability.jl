# Copyright (c) 2021, R. Dong, C. Goodbreak, H. Harrington, G. Pogudin
# Copyright (c) 2020, A. Ovchinnikov, A. Pillay, G. Pogudin, T. Scanlon

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# ------------------------------------------------------------------------------

"""
    differentiate_solution(ode, params, ic, inputs, prec)

Input:
- the same as for `power_series_solutions`

Output:
- a tuple consisting of the power series solution and a dictionary of the form `(u, v) => power series`, where `u` is a state variable
  `v` is a state or parameter, and the power series is the partial derivative of
  the function `u` w.r.t. `v` evaluated at the solution
"""
function differentiate_solution(
    ode::ODE{P},
    params::Dict{P, T},
    ic::Dict{P, T},
    inputs::Dict{P, Array{T, 1}},
    prec::Int,
) where {T <: Generic.FieldElem, P <: MPolyRingElem{T}}
    @debug "Computing the power series solution of the system"
    ps_sol = power_series_solution(ode, params, ic, inputs, prec)
    ps_ring = parent(first(values(ps_sol)))
    for p in ode.parameters
        ps_sol[p] = ps_ring(params[p])
    end

    @debug "Building the variational system at the solution"
    # Y' = AY + B
    vars = vcat(ode.x_vars, ode.parameters)
    SA = AbstractAlgebra.matrix_space(ps_ring, length(ode.x_vars), length(ode.x_vars))
    A = SA([
        eval_at_dict(derivative(ode.x_equations[vars[i]], vars[j]), ps_sol) for
        i in 1:length(ode.x_vars), j in 1:length(ode.x_vars)
    ])
    SB = AbstractAlgebra.matrix_space(ps_ring, length(ode.x_vars), length(vars))
    B = zero(SB)
    for i in 1:length(ode.x_vars)
        for j in (length(ode.x_vars) + 1):length(vars)
            B[i, j] = eval_at_dict(derivative(ode.x_equations[vars[i]], vars[j]), ps_sol)
        end
    end
    # TODO: make use of one() function (problems modulo prime)
    initial_condition =
        zero(Nemo.matrix_space(base_ring(ode.poly_ring), length(ode.x_vars), length(vars)))
    for i in 1:length(ode.x_vars)
        initial_condition[i, i] = 1
    end

    @debug "Solving the variational system and forming the output"
    sol_var_system = ps_matrix_linear_de(A, B, initial_condition, prec)
    return (
        ps_sol,
        Dict(
            (vars[i], vars[j]) => sol_var_system[i, j] for
            i in 1:length(ode.x_vars), j in 1:length(vars)
        ),
    )
end

# ------------------------------------------------------------------------------

"""
    differentiate_output(ode, params, ic, inputs, prec)

Similar to `differentiate_solution` but computes partial derivatives of prescribed outputs
returns a dictionary of the form `y_function => Dict(var => dy/dvar)` where `dy/dvar` is the derivative
of `y_function` with respect to `var`.
"""
function differentiate_output(
    ode::ODE{P},
    params::Dict{P, T},
    ic::Dict{P, T},
    inputs::Dict{P, Array{T, 1}},
    prec::Int,
) where {T <: Generic.FieldElem, P <: MPolyRingElem{T}}
    @debug "Computing partial derivatives of the solution"
    ps_sol, sol_diff = differentiate_solution(ode, params, ic, inputs, prec)
    ps_ring = parent(first(values(ps_sol)))
    for p in ode.parameters
        ps_sol[p] = ps_ring(params[p])
    end

    @debug "Evaluating the partial derivatives of the outputs"
    result = Dict()
    for (y, g) in ode.y_equations
        result[y] = Dict()
        for x in ode.x_vars
            result[y][x] = sum(
                eval_at_dict(derivative(g, xx), ps_sol) * sol_diff[(xx, x)] for
                xx in ode.x_vars
            )
        end
        for p in ode.parameters
            result[y][p] = sum(
                eval_at_dict(derivative(g, xx), ps_sol) * sol_diff[(xx, p)] for
                xx in ode.x_vars
            )
            result[y][p] += eval_at_dict(derivative(g, p), ps_sol)
        end
    end

    return result
end

# ------------------------------------------------------------------------------

"""
    get_degree_and_coeffsize(f)

for `f` being a polynomial/rational function over rationals (`QQ`) returns a tuple
`(degree, max_coef_size)`
"""
function get_degree_and_coeffsize(f::MPolyRingElem{Nemo.QQFieldElem})
    if length(f) == 0
        return (0, 1)
    end
    max_coef = 1
    for c in coefficients(f)
        max_coef = max(max_coef, 2 * height_bits(c))
    end
    return (total_degree(f), max_coef)
end

function get_degree_and_coeffsize(
    f::Generic.FracFieldElem{<:MPolyRingElem{Nemo.QQFieldElem}},
)
    num_deg, num_coef = get_degree_and_coeffsize(numerator(f))
    den_deg, den_coef = get_degree_and_coeffsize(denominator(f))
    return (max(num_deg, den_deg), max(num_coef, den_coef))
end

# ------------------------------------------------------------------------------

"""
    assess_local_identifiability(ode::ODE{P}; funcs_to_check::Array{<: Any, 1}, prob_threshold::Float64=0.99, type=:SE, loglevel=Logging.Info) where P <: MPolyRingElem{Nemo.QQFieldElem}

Checks the local identifiability/observability of the functions in `funcs_to_check`. The result is correct with probability at least `prob_threshold`.

Call this function if you have a specific collection of parameters of which you would like to check local identifiability.

`type` can be either `:SE` (single-experiment identifiability) or `:ME` (multi-experiment identifiability).
If the type is `:ME`, states are not allowed to appear in the `funcs_to_check`.
"""
function assess_local_identifiability(
    ode::ODE{P};
    funcs_to_check::Array{<:Any, 1} = Array{Any, 1}(),
    prob_threshold::Float64 = 0.99,
    type = :SE,
    trbasis = nothing,
    loglevel = Logging.Info,
) where {P <: MPolyRingElem{Nemo.QQFieldElem}}
    restart_logging(loglevel = loglevel)
    reset_timings()
    with_logger(_si_logger[]) do
        return _assess_local_identifiability(
            ode,
            funcs_to_check = funcs_to_check,
            prob_threshold = prob_threshold,
            type = type,
            trbasis = trbasis,
        )
    end
end

function _assess_local_identifiability(
    ode::ODE{P};
    funcs_to_check::Array{<:Any, 1} = Array{Any, 1}(),
    prob_threshold::Float64 = 0.99,
    type = :SE,
    trbasis = nothing,
    known_ic::Array{<:Any, 1} = Array{Any, 1}(),
) where {P <: MPolyRingElem{Nemo.QQFieldElem}}
    if isempty(funcs_to_check)
        funcs_to_check = ode.parameters
        if type == :SE
            funcs_to_check = vcat(ode.x_vars, ode.parameters)
        end
    end

    # Checking whether the states appear in the ME case
    if type == :ME
        for f in funcs_to_check
            num, den = unpack_fraction(f)
            for v in vcat(vars(num), vars(den))
                if !(v in ode.parameters)
                    @error "Multi-experiment identifiability is not properly defined for the states"
                    throw(ArgumentError("State variable $v appears in $f"))
                end
            end
        end
    end

    # Computing the prime using Proposition 3.3 from https://doi.org/10.1006/jsco.2002.0532
    @debug "Computing the prime number"
    d, h = 1, 1
    for f in vcat(collect(values(ode.x_equations)), collect(values(ode.y_equations)))
        df, hf = get_degree_and_coeffsize(f)
        d = max(d, df)
        h = max(h, hf)
    end
    p_per_func = 1 - (1 - prob_threshold) / length(funcs_to_check)
    mu = ceil(1 / (1 - sqrt(p_per_func)))

    n = length(ode.x_vars)
    # if type == ME, we take into account the largest possible replica of the system to be considered
    if type == :ME
        n *= length(ode.parameters)
    end

    m = length(ode.y_vars)
    r = length(ode.u_vars)
    ell = length(ode.parameters)
    D = 4 * (n + ell)^2 * (n + m) * d
    Dprime =
        D * (2 * log(n + ell + r + 1) + log(mu * D)) +
        4 * (n + ell)^2 * ((n + m) * h + log(2 * n * D))
    Dprime = max(Dprime, 1.0)
    prime = Primes.nextprime(Int(ceil(2 * mu * Dprime)))
    @debug "The prime is $prime"
    F = Nemo.Native.GF(prime)

    @debug "Extending the model"
    ode_ext =
        add_outputs(ode, Dict("loc_aux_$i" => f for (i, f) in enumerate(funcs_to_check)))

    @debug "Reducing the system modulo prime"
    ode_red = reduce_ode_mod_p(ode_ext, prime)

    @debug "Computing the observability matrix (and, if ME, the bound)"
    prec = length(ode.x_vars) + length(ode.parameters)

    # Parameter values are the same across all the replicas
    params_vals = Dict(p => F(rand(1:prime)) for p in ode_red.parameters)

    prev_defect = length(ode.parameters)
    num_exp = 0
    # rows are the "parameters": parameters and initial conditions
    # columns are "observations": derivatives of the outputs
    Jac = zero(Nemo.matrix_space(F, length(ode.parameters), 1))
    output_derivatives = undef
    # the while loop is primarily for ME-deintifiability, it is adding replicas until the rank stabilizes
    # in the SE case, it will exit right away
    while true
        ic = Dict(x => F(rand(1:prime)) for x in ode_red.x_vars)
        inputs = Dict{Nemo.fpMPolyRingElem, Array{Nemo.fpFieldElem, 1}}(
            u => [F(rand(1:prime)) for i in 1:prec] for u in ode_red.u_vars
        )

        @debug "Computing the output derivatives"
        output_derivatives = differentiate_output(ode_red, params_vals, ic, inputs, prec)

        @debug "Building the matrices"
        newJac = vcat(Jac, zero(Nemo.matrix_space(F, length(ode.x_vars), ncols(Jac))))
        newJac = hcat(
            newJac,
            zero(Nemo.matrix_space(F, nrows(newJac), prec * length(ode.y_vars))),
        )
        xs_params = vcat(ode_red.x_vars, ode_red.parameters)
        for (i, y) in enumerate(ode.y_vars)
            y_red = str_to_var(var_to_str(y), ode_red.poly_ring)
            offset = 1 + num_exp * prec * length(ode.y_vars)
            for j in 1:prec
                for (k, p) in enumerate(ode_red.parameters)
                    newJac[k, offset + (i - 1) * prec + j] =
                        coeff(output_derivatives[y_red][p], j - 1)
                end
                for (k, x) in enumerate(ode_red.x_vars)
                    newJac[end - k + 1, offset + (i - 1) * prec + j] =
                        coeff(output_derivatives[y_red][x], j - 1)
                end
            end
        end
        if type == :SE
            Jac = newJac
            break
        end

        new_defect =
            length(ode.parameters) +
            LinearAlgebra.rank(newJac[(length(ode.parameters) + 1):end, :]) -
            LinearAlgebra.rank(newJac)
        if new_defect == prev_defect
            break
        end
        prev_defect = new_defect
        num_exp += 1
        Jac = newJac
    end

    if !isnothing(trbasis)
        @debug "Transcendence basis computation requested"
        reverted_Jac = zero(Nemo.matrix_space(F, size(Jac)[2], size(Jac)[1]))
        for i in 1:size(Jac)[1]
            for j in 1:size(Jac)[2]
                reverted_Jac[j, i] = Jac[size(Jac)[1] - i + 1, j]
            end
        end
        Nemo.rref!(reverted_Jac)
        # finding non-pivots
        _, nonpivots = select_pivots(reverted_Jac)

        # selecting the trbasis of polynomials
        trbasis_indices_param = [
            size(Jac)[1] - i + 1 for
            i in nonpivots if i > size(Jac)[1] - length(ode.parameters)
        ]
        for i in trbasis_indices_param
            push!(trbasis, ode.parameters[i])
        end

        # NB: states are in the inverted order (to be refactored...)
        trbasis_indices_states =
            [i for i in nonpivots if i <= size(Jac)[1] - length(ode.parameters)]
        for i in trbasis_indices_states
            push!(trbasis, ode.x_vars[i])
        end
        @debug "Transcendence basis $trbasis"
    end

    @debug "Computing the result"
    base_rank = LinearAlgebra.rank(Jac)
    result = OrderedDict{Any, Bool}()
    for i in 1:length(funcs_to_check)
        for (k, p) in enumerate(ode_red.parameters)
            Jac[k, 1] =
                coeff(output_derivatives[str_to_var("loc_aux_$i", ode_red.poly_ring)][p], 0)
        end
        if type == :SE
            for (k, x) in enumerate(ode_red.x_vars)
                Jac[end - k + 1, 1] = coeff(
                    output_derivatives[str_to_var("loc_aux_$i", ode_red.poly_ring)][x],
                    0,
                )
            end
        end
        result[funcs_to_check[i]] = LinearAlgebra.rank(Jac) == base_rank
    end
    # NB: the Jac contains now the derivatives of the last from `funcs_to_check`

    if type == :SE
        return result
    end
    return (result, num_exp)
end

# ------------------------------------------------------------------------------
