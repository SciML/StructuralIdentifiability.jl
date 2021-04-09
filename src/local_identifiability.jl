# Copyright (c) 2021, R. Dong, C. Goodbreak, H. Harrington, G. Pogudin
# Copyright (c) 2020, A. Ovchinnikov, A. Pillay, G. Pogudin, T. Scanlon

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), 
# to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#------------------------------------------------------------------------------

"""
    differentiate_solution(ode, params, ic, inputs, prec)

Input: the same as for power_series_solutions
Output: a tuple consisting of the power series solution and 
        a dictionary of the form (u, v) => power series, where u is a state variable
        v is a state or parameter, and the power series is the partial derivative of
        the function u w.r.t. v evaluated at the solution
"""
function differentiate_solution(
        ode::ODE{P},
        params::Dict{P, T},
        ic::Dict{P, T},
        inputs::Dict{P, Array{T, 1}},
        prec::Int
    ) where {T <: Generic.FieldElem, P <: MPolyElem{T}}
    @debug "Computing the power series solution of the system"
    ps_sol = power_series_solution(ode, params, ic, inputs, prec)
    ps_ring = parent(first(values(ps_sol)))
    for p in ode.parameters
        ps_sol[p] = ps_ring(params[p])
    end

    @debug "Building the variational system at the solution"
    # Y' = AY + B
    vars = vcat(ode.x_vars, ode.parameters)
    SA = AbstractAlgebra.MatrixSpace(ps_ring, length(ode.x_vars), length(ode.x_vars))
    A = SA([
        eval_at_dict(derivative(ode.x_equations[vars[i]], vars[j]), ps_sol)
        for i in 1:length(ode.x_vars), j in 1:length(ode.x_vars)
    ])
    SB = AbstractAlgebra.MatrixSpace(ps_ring, length(ode.x_vars), length(vars))
    B = zero(SB)
    for i in 1:length(ode.x_vars)
        for j in (length(ode.x_vars) + 1):length(vars)
            B[i, j] = eval_at_dict(derivative(ode.x_equations[vars[i]], vars[j]), ps_sol)
        end
    end
    # TODO: make use of one() function (problems modulo prime)
    initial_condition = zero(Nemo.MatrixSpace(base_ring(ode.poly_ring), length(ode.x_vars), length(vars)))
    for i in 1:length(ode.x_vars)
        initial_condition[i, i] = 1
    end
    
    @debug "Solving the variational system and forming the output"
    sol_var_system = ps_matrix_linear_de(A, B, initial_condition, prec)
    return (
        ps_sol, 
        Dict(
            (vars[i], vars[j]) => sol_var_system[i, j]
            for i in 1:length(ode.x_vars), j in 1:length(vars)
        )
    )
end

#------------------------------------------------------------------------------

"""
    differentiate_output(ode, params, ic, inputs, prec)

Similar to differentiate_solution but computes partial derivatives of a prescribed outputs
returns a dict from y's to dictionalries from vars to dy / dvar
"""
function differentiate_output(
        ode::ODE{P},
        params::Dict{P, T},
        ic::Dict{P, T},
        inputs::Dict{P, Array{T, 1}},
        prec::Int
    ) where {T <: Generic.FieldElem, P <: MPolyElem{T}}
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
            result[y][x] = sum([eval_at_dict(derivative(g, xx), ps_sol) * sol_diff[(xx, x)] for xx in ode.x_vars])
        end
        for p in ode.parameters
            result[y][p] = sum([eval_at_dict(derivative(g, xx), ps_sol) * sol_diff[(xx, p)] for xx in ode.x_vars])
            result[y][p] += eval_at_dict(derivative(g, p), ps_sol)
        end
    end

    return result 
end

#------------------------------------------------------------------------------

"""
    get_degree_and_coeffsize(f)

for f being a polynomial/rational function over QQ returns a tuple
(degree, max_coef_size)
"""
function get_degree_and_coeffsize(f::MPolyElem{Nemo.fmpq})
    if length(f) == 0
        return (0, 1)
    end
    max_coef = 1
    for c in coeffs(f)
        max_coef = max(max_coef, 2 * height_bits(c))
    end
    return (total_degree(f), max_coef)
end

function get_degree_and_coeffsize(f::Generic.Frac{<: MPolyElem{Nemo.fmpq}})
    num_deg, num_coef = get_degree_and_coeffsize(numerator(f))
    den_deg, den_coef = get_degree_and_coeffsize(denominator(f))
    return (max(num_deg, den_deg), max(num_coef, den_coef))
end


#------------------------------------------------------------------------------

"""
    assess_local_identifiability(ode, funcs_to_check, p)

Checks the local identifiability/observability of the functions in funcs_to_check
The result is correct with probability at least p
"""
function assess_local_identifiability(ode::ODE{P}, funcs_to_check::Array{<: Any, 1}, p::Float64 = 0.99) where P <: MPolyElem{Nemo.fmpq}
    # Computing the prime using Proposition 3.3 from https://doi.org/10.1006/jsco.2002.0532
    @debug "Computing the prime number"
    d, h = 1, 1
    for f in vcat(collect(values(ode.x_equations)), collect(values(ode.y_equations)))
        df, hf = get_degree_and_coeffsize(f)
        d = max(d, df)
        h = max(h, hf)
    end
    p_per_func = 1 - (1 - p) / length(funcs_to_check)
    mu = ceil(1 / (1 - sqrt(p_per_func)))
    n = length(ode.x_vars)
    m = length(ode.y_vars)
    r = length(ode.u_vars)
    ell = length(ode.parameters)
    D = 4 * (n + ell)^2 * (n + m) * d
    Dprime = D * (2 * log(n + ell + r + 1) + log(mu * D)) + 4 * (n + ell)^2 * ((n + m) * h + log(2 * n * D))
    prime = Primes.nextprime(Int(ceil(2 * mu * Dprime)))
    @debug "The prime is $prime"
    prime = 2^31 - 1
    F = Nemo.GF(prime)
 
    @debug "Extending the model"
    ode_ext = add_outputs(ode, Dict("loc_aux_$i" => f for (i, f) in enumerate(funcs_to_check)))

    @debug "Reducing the system modulo prime"
    ode_red = reduce_ode_mod_p(ode_ext, prime)
    prec = length(ode.x_vars) + length(ode.parameters)
    params_vals = Dict(p => F(rand(1:prime)) for p in ode_red.parameters)
    ic = Dict(x => F(rand(1:prime)) for x in ode_red.x_vars)
    inputs = Dict{Nemo.gfp_mpoly, Array{Nemo.gfp_elem, 1}}(u => [F(rand(1:prime)) for i in 1:prec] for u in ode_red.u_vars)

    @debug "Computing the output derivatives"
    output_derivatives = differentiate_output(ode_red, params_vals, ic, inputs, prec)

    @debug "Building the matrices"
    # +1 is for the function to assess
    Jac = zero(Nemo.MatrixSpace(F, length(ode.y_vars) * prec + 1, prec))
    xs_params = vcat(ode_red.x_vars, ode_red.parameters)
    for (i, y) in enumerate(ode.y_vars)
        y_red = str_to_var(var_to_str(y), ode_red.poly_ring)
        for j in 1:prec
            for (k, x) in enumerate(xs_params)
                Jac[(i - 1) * prec + j, k] = coeff(output_derivatives[y_red][x], j - 1)
            end
        end
    end

    @debug "Computing the result"
    base_rank = LinearAlgebra.rank(Jac)
    result = Array{Bool, 1}()
    for i in 1:length(funcs_to_check)
        for (k, x) in enumerate(xs_params)
            Jac[end, k] = coeff(output_derivatives[str_to_var("loc_aux_$i", ode_red.poly_ring)][x], 0)
        end
        push!(result, LinearAlgebra.rank(Jac) == base_rank)
    end

    return result
end

#------------------------------------------------------------------------------
