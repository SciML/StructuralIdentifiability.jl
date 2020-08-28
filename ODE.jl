using Dates
using Logging
using Oscar

include("power_series_utils.jl")

struct ODE
    poly_ring
    x_vars
    u_vars
    parameters
    equations
    
    function ODE(eqs, inputs, field = QQ)
        #Initialize ODE
        #equations is a dictionary x_i => f_i(x, u, params)

        poly_ring = parent(collect(values(eqs))[1])
        x_vars = collect(keys(eqs))
        u_vars = inputs
        parameters = filter(v -> (!(v in x_vars) && !(v in u_vars)), gens(poly_ring))
        new(poly_ring, x_vars, u_vars, parameters, eqs)
    end
end

#------------------------------------------------------------------------------

function power_series_solution(ode::ODE, param_values, initial_conditions, input_values, prec)
    new_varnames = map(string, vcat(ode.x_vars, map(v -> "$(v)_dot", ode.x_vars), ode.u_vars))

    new_ring, new_vars = PolynomialRing(base_ring(ode.poly_ring), new_varnames)
    equations = Array{RingElem, 1}()
    evaluation = Dict(k => new_ring(v) for (k, v) in param_values)
    for v in vcat(ode.x_vars, ode.u_vars)
        evaluation[v] = str_to_var(string(v), new_ring)
    end
    for (v, eq) in ode.equations
        num, den = map(p -> eval_at_dict(p, evaluation), unpack_fraction(eq))
        push!(equations, den * str_to_var("$(v)_dot", new_ring) - num)
    end
    new_inputs = Dict(str_to_var(string(k), new_ring) => v for (k, v) in input_values)
    new_ic = Dict(str_to_var(string(k), new_ring) => v for (k, v) in initial_conditions)
    result = ps_ode_solution(equations, new_ic, new_inputs, prec)
    return Dict(v => result[str_to_var(string(v), new_ring)] for v in ode.x_vars)
end
