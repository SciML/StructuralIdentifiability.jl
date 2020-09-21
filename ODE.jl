using Dates
using Logging
using Oscar

include("power_series_utils.jl")

# P is the type of polynomials in the rhs of the ODE system
struct ODE{P}
    poly_ring::MPolyRing
    x_vars::Array{P, 1}
    u_vars::Array{P, 1}
    parameters::Array{P, 1}
    equations::Dict{P, <: Union{P, Generic.Frac{P}}}
    
    function ODE{P}(eqs::Dict{P, <: Union{P, Generic.Frac{P}}}, inputs::Array{P, 1}) where {P <: MPolyElem{<: FieldElem}}
        #Initialize ODE
        #equations is a dictionary x_i => f_i(x, u, params)

        num, den = unpack_fraction(collect(values(eqs))[1])
        poly_ring = parent(num)
        x_vars = collect(keys(eqs))
        u_vars = inputs
        parameters = filter(v -> (!(v in x_vars) && !(v in u_vars)), gens(poly_ring))
        new(poly_ring, x_vars, u_vars, parameters, eqs)
    end
end

#------------------------------------------------------------------------------

function power_series_solution(
        ode::ODE{P},
        param_values::Dict{P, T},
        initial_conditions::Dict{P, T},
        input_values::Dict{P, Array{T, 1}},
        prec::Int
    ) where {T <: FieldElem, P <: MPolyElem{T}}
    """
    Input:
        - ode, an ode to solve
        - param_values, initial_conditions - parameter values and initial conditions to plug in
          both are dictionaries variable => value
        - input_values - power series for the inpiuts presented as a dictionary
          variable => list of coefficients
        - prec, the precision
    Output: computes a power series solution with precision prec presented as a dictionary
            variable => corresponding coordiante of the solution
    """
    new_varnames = map(string, vcat(ode.x_vars, map(v -> "$(v)_dot", ode.x_vars), ode.u_vars))

    new_ring, new_vars = PolynomialRing(base_ring(ode.poly_ring), new_varnames)
    equations = Array{P, 1}()
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
    return Dict(v => result[str_to_var(string(v), new_ring)] for v in vcat(ode.x_vars, ode.u_vars))
end

#------------------------------------------------------------------------------

function power_series_solution(
        ode::ODE{P},
        param_values::Dict{P, Int},
        initial_conditions::Dict{P, Int},
        input_values::Dict{P, Array{Int, 1}},
        prec::Int
    ) where P <: MPolyElem{<: FieldElem}
    bring = base_ring(ode.poly_ring)
    return power_series_solution(
        ode,
        Dict(p => bring(v) for (p, v) in param_values),
        Dict(x => bring(v) for (x, v) in initial_conditions),
        Dict(u => map(v -> bring(v), vv) for (u, vv) in input_values),
        prec
    )
end
#------------------------------------------------------------------------------

function _reduce_poly_mod_p(poly::MPolyElem{Nemo.fmpq}, p::Int)
    """
    Reduces a polynomial over Q modulo p
    """
    den = denominator(poly)
    num = change_base_ring(ZZ, den * poly)
    if GF(p)(den) == 0
        throw(Base.ArgumentError("Prime $p divides the denominator of $poly"))
    end
    return change_base_ring(GF(p), num) * (1 // GF(p)(den))
end

#--------------------------------------

function reduce_ode_mod_p(ode::ODE{<: MPolyElem{Nemo.fmpq}}, p::Int)
    """
    Input: ode is an ODE over QQ, p is a prime number
    Output: the reduction mod p, throws an exception if p divides one of the denominators
    """
    new_ring, new_vars = PolynomialRing(GF(p), map(string, gens(ode.poly_ring)))
    new_type = typeof(new_vars[1])
    new_inputs = map(u -> str_to_var(string(u), new_ring), ode.u_vars)
    new_eqs = Dict{new_type, Union{new_type, Generic.Frac{new_type}}}()
    for (v, f) in ode.equations
        new_v = str_to_var(string(v), new_ring)
        if applicable(numerator, f)
            # if f is a rational function
            num, den = map(poly -> _reduce_poly_mod_p(poly, p), [numerator(f), denominator(f)])
            if den == 0
                throw(Base.ArgumentError("Prime $p divides the denominator of $poly"))
            end
            new_eqs[new_v] = num // den
        else
            new_eqs[new_v] = _reduce_poly_mod_p(f, p)
        end
    end
    return ODE{new_type}(new_eqs, new_inputs)
end

#------------------------------------------------------------------------------

function print_for_SIAN(ode::ODE{P}, outputs::Array{P, 1}) where P <: MPolyElem{<: FieldElem}
    """
    Prints the ODE in the format accepted by SIAN (https://github.com/pogudingleb/SIAN)
    """
    var_to_str = Dict(x => "$(x)(t)" for x in vcat(ode.x_vars, ode.u_vars))
    merge!(var_to_str, Dict(p => string(p) for p in ode.parameters))
    R_print, vars_print = PolynomialRing(base_ring(ode.poly_ring), [var_to_str[v] for v in gens(ode.poly_ring)])
    result = ""

    function _lhs_to_str(lhs)
        num, den = unpack_fraction(lhs)
        result = string(evaluate(num, vars_print))
        if den != 1
             result = "($result) / ($(evaluate(den, vars_print)))"
        end
        return result
    end

    for (x, f) in ode.equations
        result = result * "diff($(x)(t), t) = $(_lhs_to_str(f)), \n"
    end
    for (y_ind, g) in enumerate(outputs)
        result = result * "y_var_$y_ind(t) = $(_lhs_to_str(g)), \n"
    end
    return result
end

#------------------------------------------------------------------------------
