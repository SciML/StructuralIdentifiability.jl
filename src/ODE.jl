using Dates
using Logging
using MacroTools
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
    new_varnames = map(var_to_str, vcat(ode.x_vars, ode.u_vars))
    append!(new_varnames, map(v -> var_to_str(v) * "_dot", ode.x_vars))

    new_ring, new_vars = PolynomialRing(base_ring(ode.poly_ring), new_varnames)
    equations = Array{P, 1}()
    evaluation = Dict(k => new_ring(v) for (k, v) in param_values)
    for v in vcat(ode.x_vars, ode.u_vars)
        evaluation[v] = switch_ring(v, new_ring)
    end
    for (v, eq) in ode.equations
        num, den = map(p -> eval_at_dict(p, evaluation), unpack_fraction(eq))
        push!(equations, den * str_to_var(var_to_str(v) * "_dot", new_ring) - num)
    end
    new_inputs = Dict(switch_ring(k, new_ring) => v for (k, v) in input_values)
    new_ic = Dict(switch_ring(k, new_ring) => v for (k, v) in initial_conditions)
    result = ps_ode_solution(equations, new_ic, new_inputs, prec)
    return Dict(v => result[switch_ring(v, new_ring)] for v in vcat(ode.x_vars, ode.u_vars))
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
    new_ring, new_vars = PolynomialRing(GF(p), map(var_to_str, gens(ode.poly_ring)))
    new_type = typeof(new_vars[1])
    new_inputs = map(u -> switch_ring(u, new_ring), ode.u_vars)
    new_eqs = Dict{new_type, Union{new_type, Generic.Frac{new_type}}}()
    for (v, f) in ode.equations
        new_v = switch_ring(v, new_ring)
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
    var_to_str = Dict(x => var_to_str(x) * "(t)" for x in vcat(ode.x_vars, ode.u_vars))
    merge!(var_to_str, Dict(p => var_to_str(p) for p in ode.parameters))
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
        result = result * "diff(" * var_to_str(x) * "(t), t) = $(_lhs_to_str(f)), \n"
    end
    for (y_ind, g) in enumerate(outputs)
        result = result * "y_var_$y_ind(t) = $(_lhs_to_str(g)), \n"
    end
    return result
end

#------------------------------------------------------------------------------

function macrohelper_extract_vars(equations::Array{Expr, 1}, summary::Bool=true)
    funcs, x_vars, all_symb = Set(), Set(), Set()
    aux_symb = Set([:(+), :(-), :(=), :(*), :(^), :t, :(/), :(//)])
    for eq in equations
        MacroTools.postwalk(
            x -> begin 
                if @capture(x, f_'(t)) 
                    push!(x_vars, f)
                    push!(all_symb, f)
                elseif @capture(x, f_(t))
                    push!(funcs, f)
                elseif (x isa Symbol) && !(x in aux_symb)
                    push!(all_symb, x)
                end
                return x
            end, 
            eq
        )
    end
    u_vars = setdiff(funcs, x_vars)
    params = setdiff(all_symb, union(x_vars, u_vars))
    all_symb = collect(all_symb)
    if summary
        print("Summary of the model:\n")
        print("State variables: ", join(map(string, collect(x_vars)), ", "), "\n")
        print("Parameter: ", join(map(string, collect(params)), ", "), "\n")
        print("Inputs: ", join(map(string, collect(u_vars)), ", "), "\n")
    end
    return collect(u_vars), collect(all_symb)
end

#------------------------------------------------------------------------------

function macrohelper_clean(ex::Expr)
    ex = MacroTools.postwalk(x -> @capture(x, f_'(t)) ? f : x, ex)
    ex = MacroTools.postwalk(x -> @capture(x, f_(t)) ? f : x, ex)
    ex = MacroTools.postwalk(x -> x == :(/) ? :(//) : x, ex)
    return ex
end

#------------------------------------------------------------------------------

macro ODEmodel(ex::Expr...)
    """
    Macros for creating an ODE from a list of equations
    Also injects all variables into the global scope
    """
    extra_params = []
    equations = [ex...]
    if equations[end].head == :vect
        extra_params = equations[end].args
        equations = equations[1:end - 1]
    end

    u_vars, all_symb = macrohelper_extract_vars(equations)
    append!(all_symb, extra_params)
    
    # creating the polynomial ring
    vars_list = :([$(all_symb...)])
    R = gensym()
    vars_aux = gensym()
    exp_ring = :(($R, $vars_aux) = PolynomialRing(QQ, map(string, $all_symb)))
    assignments = [:($(all_symb[i]) = $vars_aux[$i]) for i in 1:length(all_symb)]
    
    # preparing equations
    equations = map(macrohelper_clean, equations)
    dict_var = gensym()
    dict_create_expr = :($dict_var = Dict{fmpq_mpoly, Union{fmpq_mpoly, Generic.Frac{fmpq_mpoly}}}())
    eqs_expr = []
    for eq in equations
        if eq.head != :(=)
            throw("Problem with parsing at $eq") 
        end
        lhs, rhs = eq.args[1:2]
        loc_all_symb = macrohelper_extract_vars([rhs], false)[2]
        if isempty(loc_all_symb)
            push!(eqs_expr, :($dict_var[$lhs] = $R($rhs)))
        else
            push!(eqs_expr, :($dict_var[$lhs] = ($rhs)))
        end
    end
    
    # creating the ode object
    ode_expr = :(ODE{fmpq_mpoly}($dict_var, Array{fmpq_mpoly}([$(u_vars...)])))
    
    result = Expr(
        :block, 
        exp_ring, assignments..., 
        dict_create_expr, eqs_expr..., 
        ode_expr
    )
    return esc(result)
end

#------------------------------------------------------------------------------
