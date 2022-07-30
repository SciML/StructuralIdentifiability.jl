# P is the type of polynomials in the rhs of the ODE system
"""
The main structure that represents input ODE system.

Stores information about states (`x_vars`), outputs (`y_vars`), inputs (`u_vars`), parameters (`parameters`) and the equations.

This structure is constructed via `@ODEmodel` macro.
"""
struct ODE{P}
    poly_ring::MPolyRing
    x_vars::Array{P, 1}
    y_vars::Array{P, 1}
    u_vars::Array{P, 1}
    parameters::Array{P, 1}
    x_equations::Dict{P, <: Union{P, Generic.Frac{P}}}
    y_equations::Dict{P, <: Union{P, Generic.Frac{P}}}
    function ODE{P}(
            x_eqs::Dict{P, <: Union{P, Generic.Frac{P}}}, 
            y_eqs::Dict{P, <: Union{P, Generic.Frac{P}}},    
            inputs::Array{P, 1}
        ) where {P <: MPolyElem{<: FieldElem}}
        # Initialize ODE
        # x_eqs is a dictionary x_i => f_i(x, u, params)
        # y_eqs is a dictionary y_i => g_i(x, u, params)

        num, den = unpack_fraction(collect(values(x_eqs))[1])
        poly_ring = parent(num)
        if !all(isascii.(string.(gens(poly_ring))))
            nonascii_chars = filter(g->!isascii(g), string.(gens(poly_ring)))
            st = join(nonascii_chars, ", ")
            @warn "Non-ascii characters are not supported by Singular: " * st
        end
        x_vars = collect(keys(x_eqs))
        y_vars = collect(keys(y_eqs))
        u_vars = inputs
        parameters = filter(v -> (!(v in x_vars) && !(v in u_vars) && !(v in y_vars)), gens(poly_ring))
        new{P}(poly_ring, x_vars, y_vars, u_vars, parameters, x_eqs, y_eqs)
    end
end

#------------------------------------------------------------------------------

function add_outputs(ode::ODE{P}, extra_y::Dict{String, <: RingElem}) where P <: MPolyElem
    new_var_names = vcat(collect(map(var_to_str, gens(ode.poly_ring))), collect(keys(extra_y)))
    new_ring, new_vars = Nemo.PolynomialRing(base_ring(ode.poly_ring), new_var_names)
    
    new_x_eqs = Dict{P, Union{P, Generic.Frac{P}}}(parent_ring_change(x, new_ring) => parent_ring_change(f, new_ring) for (x, f) in ode.x_equations)
    new_y_eqs = Dict{P, Union{P, Generic.Frac{P}}}(parent_ring_change(y, new_ring) => parent_ring_change(g, new_ring) for (y, g) in ode.y_equations)
    extra_y_eqs = Dict{P, Union{P, Generic.Frac{P}}}(str_to_var(y, new_ring) => parent_ring_change(g, new_ring) for (y, g) in extra_y)
    merge!(new_y_eqs, extra_y_eqs)
    new_us = map(v -> switch_ring(v, new_ring), ode.u_vars)
    return ODE{P}(new_x_eqs, new_y_eqs, new_us)
end

#------------------------------------------------------------------------------

"""
    set_parameter_values(ode, param_values)

Input:
- `ode` - an ODE as above
- `param_values` - values for (possibly, some of) the parameters as dictionary `parameter` => `value`

Output: 
- new ode with the parameters in param_values plugged with the given numbers
"""
function set_parameter_values(ode::ODE{P}, param_values::Dict{P, T}) where {T <: FieldElem, P <: MPolyElem{T}}
    new_vars = map(var_to_str, [v for v in gens(ode.poly_ring) if !(v in keys(param_values))])
    small_ring, small_vars = Nemo.PolynomialRing(base_ring(ode.poly_ring), new_vars)
    eval_dict = Dict(str_to_var(v, ode.poly_ring) => str_to_var(v, small_ring) for v in new_vars)
    merge!(eval_dict, Dict(p => small_ring(val) for (p, val) in param_values))

    return ODE{P}(
        Dict{P, Union{P, Generic.Frac{P}}}(eval_at_dict(v, eval_dict) => eval_at_dict(f, eval_dict) for (v, f) in ode.x_equations),
        Dict{P, Union{P, Generic.Frac{P}}}(eval_at_dict(v, eval_dict) => eval_at_dict(f, eval_dict) for (v, f) in ode.y_equations),
        [eval_at_dict(u, eval_dict) for u in ode.u_vars]
    )
end

#------------------------------------------------------------------------------

"""
    power_series_solution(ode, param_values, initial_conditions, input_values, prec)

Input:
- `ode` - an ode to solve
- `param_values` - parameter values, must be a dictionary mapping parameter to a value
- `initial_conditions` - initial conditions of `ode`, must be a dictionary mapping state variable to a value
- `input_values` - power series for the inpiuts presented as a dictionary variable => list of coefficients
- `prec` - the precision of solutions

Output: 
- computes a power series solution with precision prec presented as a dictionary variable => corresponding coordiante of the solution
"""
function power_series_solution(
        ode::ODE{P},
        param_values::Dict{P, T},
        initial_conditions::Dict{P, T},
        input_values::Dict{P, Array{T, 1}},
        prec::Int
    ) where {T <: FieldElem, P <: MPolyElem{T}}
    new_varnames = map(var_to_str, vcat(ode.x_vars, ode.u_vars))
    append!(new_varnames, map(v -> var_to_str(v) * "_dot", ode.x_vars))

    new_ring, new_vars = Nemo.PolynomialRing(base_ring(ode.poly_ring), new_varnames)
    equations = Array{P, 1}()
    evaluation = Dict(k => new_ring(v) for (k, v) in param_values)
    for v in vcat(ode.x_vars, ode.u_vars)
        evaluation[v] = switch_ring(v, new_ring)
    end
    for (v, eq) in ode.x_equations
        num, den = map(p -> eval_at_dict(p, evaluation), unpack_fraction(eq))
        push!(equations, den * str_to_var(var_to_str(v) * "_dot", new_ring) - num)
    end
    new_inputs = Dict(switch_ring(k, new_ring) => v for (k, v) in input_values)
    new_ic = Dict(switch_ring(k, new_ring) => v for (k, v) in initial_conditions)
    result = ps_ode_solution(equations, new_ic, new_inputs, prec)

    # Evaluation outputs
    result = Dict(v => result[switch_ring(v, new_ring)] for v in vcat(ode.x_vars, ode.u_vars))
    ps_ring = parent(first(values(result)))
    for p in ode.parameters
        result[p] = ps_ring(param_values[p])
    end
    eval_outputs = []
    for (y, g) in ode.y_equations
        result[y] = eval_at_dict(g, result)
    end
    return result
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

"""
    _reduce_mod_p(f, p)

Reduces a polynomial/rational function over Q modulo p
"""
function _reduce_mod_p(poly::fmpq_mpoly, p::Int)
    den = denominator(poly)
    num = change_base_ring(Nemo.ZZ, den * poly)
    if Nemo.GF(p)(den) == 0
        throw(Base.ArgumentError("Prime $p divides the denominator of $poly"))
    end
    return change_base_ring(Nemo.GF(p), num) * (1 // Nemo.GF(p)(den))
end

function _reduce_mod_p(rat::Generic.Frac{fmpq_mpoly}, p::Int)
    num, den = map(poly -> _reduce_mod_p(poly, p), [numerator(rat), denominator(rat)])
    if den == 0
        throw(Base.ArgumentError("Prime $p divides the denominator of $rat"))
    end
    return num // den
end

#--------------------------------------


"""
    reduce_ode_mod_p(ode, p)

Input: ode is an ODE over QQ, p is a prime number
Output: the reduction mod p, throws an exception if p divides one of the denominators
"""
function reduce_ode_mod_p(ode::ODE{<: MPolyElem{Nemo.fmpq}}, p::Int)
    new_ring, new_vars = Nemo.PolynomialRing(Nemo.GF(p), map(var_to_str, gens(ode.poly_ring)))
    new_type = typeof(new_vars[1])
    new_inputs = map(u -> switch_ring(u, new_ring), ode.u_vars)
    new_x_eqs = Dict{new_type, Union{new_type, Generic.Frac{new_type}}}()
    new_y_eqs = Dict{new_type, Union{new_type, Generic.Frac{new_type}}}()
    for (old, new) in Dict(ode.x_equations => new_x_eqs, ode.y_equations => new_y_eqs)
        for (v, f) in old
            new_v = switch_ring(v, new_ring)
            new[new_v] = _reduce_mod_p(f, p)
        end
    end
    return ODE{new_type}(new_x_eqs, new_y_eqs, new_inputs)
end

#------------------------------------------------------------------------------

function macrohelper_extract_vars(equations::Array{Expr, 1}, ders_ok::Bool=false)
    funcs, x_vars, all_symb = Set(), Set(), Set()
    aux_symb = Set([:(+), :(-), :(=), :(*), :(^), :t, :(/), :(//)])
    for eq in equations
        MacroTools.postwalk(
            x -> begin 
                if @capture(x, f_'(t))
		    if !ders_ok
                        throw(Base.ArgumentError("Derivative are not allowed in the right-hand side"))
		    end
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
    io_vars = setdiff(funcs, x_vars)
    all_symb = collect(all_symb)
    return collect(x_vars), collect(io_vars), collect(all_symb)
end

function macrohelper_extract_vars(equations::Array{Symbol, 1})
    return macrohelper_extract_vars(map(Expr, equations))
end

#------------------------------------------------------------------------------

function macrohelper_clean(ex::Expr)
    ex = MacroTools.postwalk(x -> @capture(x, f_'(t)) ? f : x, ex)
    ex = MacroTools.postwalk(x -> @capture(x, f_(t)) ? f : x, ex)
    ex = MacroTools.postwalk(x -> x == :(/) ? :(//) : x, ex)
    return ex
end

#------------------------------------------------------------------------------

""" 
Macro for creating an ODE from a list of equations.
Also injects all variables into the global scope.

This macro accepts a sybolically written ODE system and generates an `ODE` structure instance:
```julia
ode = @ODEmodel(
    x1'(t) = -k1 * x1(t),
    y1(t) = x1(t)
)
```
"""
macro ODEmodel(ex::Expr...)
    equations = [ex...]
    x_vars, io_vars, all_symb = macrohelper_extract_vars(equations; ders_ok=True)
   
    # creating the polynomial ring
    vars_list = :([$(all_symb...)])
    R = gensym()
    vars_aux = gensym()
    exp_ring = :(($R, $vars_aux) = StructuralIdentifiability.Nemo.PolynomialRing(
        StructuralIdentifiability.Nemo.QQ, 
        map(string, $all_symb)
    ))
    assignments = [:($(all_symb[i]) = $vars_aux[$i]) for i in 1:length(all_symb)]
    
    # preparing equations
    equations = map(macrohelper_clean, equations)
    x_dict = gensym()
    y_dict = gensym()
    y_vars = Set()
    x_dict_create_expr = :($x_dict = Dict{
        StructuralIdentifiability.Nemo.fmpq_mpoly, 
        Union{StructuralIdentifiability.Nemo.fmpq_mpoly, StructuralIdentifiability.AbstractAlgebra.Generic.Frac{StructuralIdentifiability.Nemo.fmpq_mpoly}}
    }())
    y_dict_create_expr = :($y_dict = Dict{
        StructuralIdentifiability.Nemo.fmpq_mpoly, 
        Union{StructuralIdentifiability.Nemo.fmpq_mpoly, StructuralIdentifiability.AbstractAlgebra.Generic.Frac{StructuralIdentifiability.Nemo.fmpq_mpoly}}
    }())
    eqs_expr = []
    for eq in equations
        if eq.head != :(=)
            throw("Problem with parsing at $eq") 
        end
        lhs, rhs = eq.args[1:2]
        loc_all_symb = macrohelper_extract_vars([rhs])[3]
        to_insert = undef
        if lhs in x_vars
            to_insert = x_dict
        elseif lhs in io_vars
            to_insert = y_dict
            push!(y_vars, lhs)
        else
            throw("Unknown left-hand side $lhs")
        end
        if isempty(loc_all_symb)
            push!(eqs_expr, :($to_insert[$lhs] = $R($rhs)))
        else
            push!(eqs_expr, :($to_insert[$lhs] = ($rhs)))
        end
    end

    u_vars = setdiff(io_vars, y_vars)
    params = setdiff(all_symb, union(x_vars, y_vars, u_vars))
    @info "Summary of the model:"
    @info "State variables: " * join(map(string, collect(x_vars)), ", ")
    @info "Parameters: " * join(map(string, collect(params)), ", ")
    @info "Inputs: " * join(map(string, collect(u_vars)), ", ")
    @info "Outputs: " * join(map(string, collect(y_vars)), ", ")
   
    # creating the ode object
    ode_expr = :(ODE{StructuralIdentifiability.Nemo.fmpq_mpoly}($x_dict, $y_dict, Array{StructuralIdentifiability.Nemo.fmpq_mpoly}([$(u_vars...)])))
    
    result = Expr(
        :block,
        exp_ring, assignments..., 
        x_dict_create_expr, y_dict_create_expr, eqs_expr..., 
        ode_expr
    )
    return esc(result)
end

#------------------------------------------------------------------------------

function Base.show(io::IO, ode::ODE)
    varstr = Dict(x => var_to_str(x) * "(t)" for x in vcat(ode.x_vars, ode.u_vars, ode.y_vars))
    merge!(varstr, Dict(p => var_to_str(p) for p in ode.parameters))
    R_print, vars_print = Nemo.PolynomialRing(base_ring(ode.poly_ring), [varstr[v] for v in gens(ode.poly_ring)])
    for (x, eq) in ode.x_equations
        print(io, var_to_str(x) * "'(t) = ")
        print(io, evaluate(eq, vars_print))
        print(io, "\n")
    end
    for (y, eq) in ode.y_equations
        print(io, var_to_str(y) * "(t) = ")
        print(io, evaluate(eq, vars_print))
        print(io, "\n")
    end
end

#------------------------------------------------------------------------------
"""
    function PreprocessODE(de::ModelingToolkit.ODESystem, measured_quantities::Array{ModelingToolkit.Equation})
    
Input:
- `de` - ModelingToolkit.ODESystem, a system for identifiability query
- `measured_quantities` - array of output functions

Output: 
- `ODE` object containing required data for identifiability assessment
"""
function PreprocessODE(de::ModelingToolkit.ODESystem, measured_quantities::Array{ModelingToolkit.Equation})
    @info "Preproccessing `ModelingToolkit.ODESystem` object"
    diff_eqs = filter(eq->!(ModelingToolkit.isoutput(eq.lhs)), ModelingToolkit.equations(de))
    y_functions = [each.lhs for each in measured_quantities]
    inputs = filter(v->ModelingToolkit.isinput(v), ModelingToolkit.states(de))
    state_vars = filter(s->!(ModelingToolkit.isinput(s) || ModelingToolkit.isoutput(s)), ModelingToolkit.states(de))
    params = ModelingToolkit.parameters(de)
    t = ModelingToolkit.arguments(measured_quantities[1].lhs)[1]
    params_from_measured_quantities = ModelingToolkit.parameters(ModelingToolkit.ODESystem(measured_quantities, t, name=:DataSeries))
    params = union(params, params_from_measured_quantities)
    
    input_symbols = vcat(state_vars, y_functions, inputs, params)
    generators = string.(input_symbols)
    generators = map(g->replace(g, "(t)"=>""), generators)
    R, gens_ = Nemo.PolynomialRing(Nemo.QQ, generators)
    state_eqn_dict = Dict{StructuralIdentifiability.Nemo.fmpq_mpoly,Union{StructuralIdentifiability.Nemo.fmpq_mpoly,StructuralIdentifiability.Nemo.Generic.Frac{StructuralIdentifiability.Nemo.fmpq_mpoly}}}()
    out_eqn_dict = Dict{StructuralIdentifiability.Nemo.fmpq_mpoly,Union{StructuralIdentifiability.Nemo.fmpq_mpoly,StructuralIdentifiability.Nemo.Generic.Frac{StructuralIdentifiability.Nemo.fmpq_mpoly}}}()
    
    for i in 1:length(diff_eqs)
        if !(typeof(diff_eqs[i].rhs) <: Number)
            state_eqn_dict[substitute(state_vars[i], input_symbols.=>gens_)] = eval_at_nemo(diff_eqs[i].rhs, Dict(input_symbols.=>gens_))
        else
            state_eqn_dict[substitute(state_vars[i], input_symbols.=>gens_)] = R(diff_eqs[i].rhs) 
        end
    end
    for i in 1:length(measured_quantities)
        out_eqn_dict[substitute(y_functions[i], input_symbols.=> gens_)] = eval_at_nemo(measured_quantities[i].rhs, Dict(input_symbols.=>gens_))
    end
    
    inputs_ = [substitute(each,  input_symbols .=> gens_) for each in inputs]
    if isequal(length(inputs_), 0)
        inputs_ = Vector{StructuralIdentifiability.Nemo.fmpq_mpoly}()
    end
    return (StructuralIdentifiability.ODE{StructuralIdentifiability.Nemo.fmpq_mpoly}(state_eqn_dict, out_eqn_dict, inputs_), input_symbols, gens_)
end
