"""
The main structure that represents input ODE system.

Stores information about states (`x_vars`), outputs (`y_vars`), inputs (`u_vars`), parameters (`parameters`) and the equations.

This structure is constructed via `@ODEmodel` macro.
"""
struct ODE{P} # P is the type of polynomials in the rhs of the ODE system
    poly_ring::MPolyRing
    x_vars::Array{P, 1}
    y_vars::Array{P, 1}
    u_vars::Array{P, 1}
    parameters::Array{P, 1}
    x_equations::Dict{P, <:Union{P, Generic.Frac{P}}}
    y_equations::Dict{P, <:Union{P, Generic.Frac{P}}}

    function ODE{P}(
        x_vars::Array{P, 1},
        y_vars::Array{P, 1},
        x_eqs::Dict{P, <:Union{P, Generic.Frac{P}}},
        y_eqs::Dict{P, <:Union{P, Generic.Frac{P}}},
        inputs::Array{P, 1},
    ) where {P <: MPolyElem{<:FieldElem}}
        # Initialize ODE
        # x_eqs is a dictionary x_i => f_i(x, u, params)
        # y_eqs is a dictionary y_i => g_i(x, u, params)
        if isempty(y_eqs)
            @info "Could not find output variables in the model."
        end
        poly_ring = parent(first(vcat(y_vars, x_vars)))
        u_vars = inputs
        parameters = filter(
            v -> (!(v in x_vars) && !(v in u_vars) && !(v in y_vars)),
            gens(poly_ring),
        )
        new{P}(poly_ring, x_vars, y_vars, u_vars, parameters, x_eqs, y_eqs)
    end

    function ODE{P}(
        x_eqs::Dict{P, <:Union{P, Generic.Frac{P}}},
        y_eqs::Dict{P, <:Union{P, Generic.Frac{P}}},
        inputs::Array{P, 1},
    ) where {P <: MPolyElem{<:FieldElem}}
        x_vars = collect(keys(x_eqs))
        y_vars = collect(keys(y_eqs))
        return ODE{P}(x_vars, y_vars, x_eqs, y_eqs, inputs)
    end
end

function Base.parent(ode::ODE)
    return ode.poly_ring
end

#------------------------------------------------------------------------------

function add_outputs(ode::ODE{P}, extra_y::Dict{String, <:RingElem}) where {P <: MPolyElem}
    new_var_names =
        vcat(collect(map(var_to_str, gens(ode.poly_ring))), collect(keys(extra_y)))
    new_ring, new_vars = Nemo.PolynomialRing(base_ring(ode.poly_ring), new_var_names)

    new_x = Array{P, 1}([parent_ring_change(x, new_ring) for x in ode.x_vars])
    new_x_eqs = Dict{P, Union{P, Generic.Frac{P}}}(
        parent_ring_change(x, new_ring) => parent_ring_change(f, new_ring) for
        (x, f) in ode.x_equations
    )
    new_y = Array{P, 1}([parent_ring_change(y, new_ring) for y in ode.y_vars])
    for y in keys(extra_y)
        push!(new_y, str_to_var(y, new_ring))
    end
    new_y_eqs = Dict{P, Union{P, Generic.Frac{P}}}(
        parent_ring_change(y, new_ring) => parent_ring_change(g, new_ring) for
        (y, g) in ode.y_equations
    )
    extra_y_eqs = Dict{P, Union{P, Generic.Frac{P}}}(
        str_to_var(y, new_ring) => parent_ring_change(g, new_ring) for (y, g) in extra_y
    )
    merge!(new_y_eqs, extra_y_eqs)
    new_us = map(v -> switch_ring(v, new_ring), ode.u_vars)
    return ODE{P}(new_x, new_y, new_x_eqs, new_y_eqs, new_us)
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
function set_parameter_values(
    ode::ODE{P},
    param_values::Dict{P, T},
) where {T <: FieldElem, P <: MPolyElem{T}}
    new_vars =
        map(var_to_str, [v for v in gens(ode.poly_ring) if !(v in keys(param_values))])
    small_ring, small_vars = Nemo.PolynomialRing(base_ring(ode.poly_ring), new_vars)
    eval_dict =
        Dict(str_to_var(v, ode.poly_ring) => str_to_var(v, small_ring) for v in new_vars)
    merge!(eval_dict, Dict(p => small_ring(val) for (p, val) in param_values))

    return ODE{P}(
        Dict{P, Union{P, Generic.Frac{P}}}(
            eval_at_dict(v, eval_dict) => eval_at_dict(f, eval_dict) for
            (v, f) in ode.x_equations
        ),
        Dict{P, Union{P, Generic.Frac{P}}}(
            eval_at_dict(v, eval_dict) => eval_at_dict(f, eval_dict) for
            (v, f) in ode.y_equations
        ),
        [eval_at_dict(u, eval_dict) for u in ode.u_vars],
    )
end

function _to_rational(x::Float64)
    s = "$x"
    point = first(findfirst(".", s))
    numer = parse(BigInt, s[1:(point - 1)] * s[(point + 1):end])
    denom = 10^(length(s) - point)
    return Nemo.QQ(numer // denom)
end

function _to_rational(x)
    return Nemo.QQ(x)
end

function set_parameter_values(
    ode::ODE{P},
    param_values::Dict{P, <:Number},
) where {P <: MPolyElem}
    new_values = Dict{P, fmpq}(x => _to_rational(v) for (x, v) in param_values)
    return set_parameter_values(ode, new_values)
end

#------------------------------------------------------------------------------

"""
    power_series_solution(ode, param_values, initial_conditions, input_values, prec)

Input:
- `ode` - an ode to solve
- `param_values` - parameter values, must be a dictionary mapping parameter to a value
- `initial_conditions` - initial conditions of `ode`, must be a dictionary mapping state variable to a value
- `input_values` - power series for the inputs presented as a dictionary variable => list of coefficients
- `prec` - the precision of solutions

Output:
- computes a power series solution with precision prec presented as a dictionary variable => corresponding coordinate of the solution
"""
function power_series_solution(
    ode::ODE{P},
    param_values::Dict{P, T},
    initial_conditions::Dict{P, T},
    input_values::Dict{P, Array{T, 1}},
    prec::Int,
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
    result =
        Dict(v => result[switch_ring(v, new_ring)] for v in vcat(ode.x_vars, ode.u_vars))
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
    prec::Int,
) where {P <: MPolyElem{<:FieldElem}}
    bring = base_ring(ode.poly_ring)
    return power_series_solution(
        ode,
        Dict(p => bring(v) for (p, v) in param_values),
        Dict(x => bring(v) for (x, v) in initial_conditions),
        Dict(u => map(v -> bring(v), vv) for (u, vv) in input_values),
        prec,
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
function reduce_ode_mod_p(ode::ODE{<:MPolyElem{Nemo.fmpq}}, p::Int)
    new_ring, new_vars =
        Nemo.PolynomialRing(Nemo.GF(p), map(var_to_str, gens(ode.poly_ring)))
    new_type = typeof(new_vars[1])
    new_inputs = map(u -> switch_ring(u, new_ring), ode.u_vars)
    new_x = map(x -> switch_ring(x, new_ring), ode.x_vars)
    new_y = map(y -> switch_ring(y, new_ring), ode.y_vars)
    new_x_eqs = Dict{new_type, Union{new_type, Generic.Frac{new_type}}}()
    new_y_eqs = Dict{new_type, Union{new_type, Generic.Frac{new_type}}}()
    for (old, new) in Dict(ode.x_equations => new_x_eqs, ode.y_equations => new_y_eqs)
        for (v, f) in old
            new_v = switch_ring(v, new_ring)
            new[new_v] = _reduce_mod_p(f, p)
        end
    end
    return ODE{new_type}(new_x, new_y, new_x_eqs, new_y_eqs, new_inputs)
end

#------------------------------------------------------------------------------

function _extract_aux!(funcs, all_symb, eq, ders_ok = false)
    aux_symb = Set([:(+), :(-), :(=), :(*), :(^), :t, :(/), :(//)])
    MacroTools.postwalk(
        x -> begin
            if @capture(x, f_'(t))
                if !ders_ok
                    throw(
                        Base.ArgumentError(
                            "Derivative are not allowed in the right-hand side",
                        ),
                    )
                end
                push!(all_symb, f)
            elseif @capture(x, f_(t))
                push!(funcs, f)
            elseif (x isa Symbol) && !(x in aux_symb)
                push!(all_symb, x)
            end
            return x
        end,
        eq,
    )
end

"""
  For an expression of the form f'(t) or f(t) returns (f, true) and (f, false), resp
"""
function _get_var(expr)
    if @capture(expr, f_'(t))
        return (f, true)
    end
    if @capture(expr, f_(t))
        return (f, false)
    end
    error("cannot extract the single function name from $expr")
end

function macrohelper_extract_vars(equations::Array{Expr, 1})
    funcs, all_symb = Set(), Set()
    x_vars, y_vars = Vector(), Vector()
    aux_symb = Set([:(+), :(-), :(=), :(*), :(^), :t, :(/), :(//)])
    for eq in equations
        if eq.head != :(=)
            _extract_aux!(funcs, all_symb, eq)
        else
            lhs, rhs = eq.args[1:2]
            _extract_aux!(funcs, all_symb, lhs, true)
            _extract_aux!(funcs, all_symb, rhs)
            (v, is_state) = _get_var(lhs)
            if is_state
                push!(x_vars, v)
            else
                push!(y_vars, v)
            end
        end
    end
    u_vars = setdiff(funcs, vcat(x_vars, y_vars))
    all_symb = collect(all_symb)
    return x_vars, y_vars, collect(u_vars), collect(all_symb)
end

function macrohelper_extract_vars(equations::Array{Symbol, 1})
    return macrohelper_extract_vars(map(Expr, equations))
end

#------------------------------------------------------------------------------

function macrohelper_clean(ex::Expr)
    ex = MacroTools.postwalk(x -> @capture(x, f_'(t)) ? f : x, ex)
    ex = MacroTools.postwalk(x -> @capture(x, f_(t)) ? f : x, ex)
    ex = MacroTools.postwalk(x -> x == :(/) ? :(//) : x, ex)
    ex = MacroTools.postwalk(x -> x isa Float64 ? rationalize(x) : x, ex)
    return ex
end

#------------------------------------------------------------------------------

"""
    macro ODEmodel

Macro for creating an ODE from a list of equations.
It also injects all variables into the global scope.

## Example

Creating a simple `ODE`:

```jldoctest
using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = a * x1(t) + u(t),
    x2'(t) = b * x2(t) + c*x1(t)*x2(t),
    y(t) = x1(t)
)
```

Here,
- `x1`, `x2` are state variables
- `y` is an output variable
- `u` is an input variable
- `a`, `b`, `c` are time-indepdendent parameters

"""
macro ODEmodel(ex::Expr...)
    equations = [ex...]
    x_vars, y_vars, u_vars, all_symb = macrohelper_extract_vars(equations)
    time_dependent = vcat(x_vars, y_vars, u_vars)
    params = sort([s for s in all_symb if !(s in time_dependent)])
    all_symb_no_t = vcat(time_dependent, params)
    all_symb_with_t = vcat([:($s(t)) for s in time_dependent], params)

    # creating the polynomial ring
    vars_list = :([$(all_symb_with_t...)])
    R = gensym()
    vars_aux = gensym()
    exp_ring = :(
        ($R, $vars_aux) = StructuralIdentifiability.Nemo.PolynomialRing(
            StructuralIdentifiability.Nemo.QQ,
            map(string, $all_symb_with_t),
        )
    )
    assignments = [:($(all_symb_no_t[i]) = $vars_aux[$i]) for i in 1:length(all_symb_no_t)]

    # setting x_vars and y_vars in the right order
    vx = gensym()
    vy = gensym()
    x_var_expr = :($vx = Vector{StructuralIdentifiability.Nemo.fmpq_mpoly}([$(x_vars...)]))
    y_var_expr = :($vy = Vector{StructuralIdentifiability.Nemo.fmpq_mpoly}([$(y_vars...)]))

    # preparing equations
    equations = map(macrohelper_clean, equations)
    x_dict = gensym()
    y_dict = gensym()
    x_dict_create_expr = :(
        $x_dict = Dict{
            StructuralIdentifiability.Nemo.fmpq_mpoly,
            Union{
                StructuralIdentifiability.Nemo.fmpq_mpoly,
                StructuralIdentifiability.AbstractAlgebra.Generic.Frac{
                    StructuralIdentifiability.Nemo.fmpq_mpoly,
                },
            },
        }()
    )
    y_dict_create_expr = :(
        $y_dict = Dict{
            StructuralIdentifiability.Nemo.fmpq_mpoly,
            Union{
                StructuralIdentifiability.Nemo.fmpq_mpoly,
                StructuralIdentifiability.AbstractAlgebra.Generic.Frac{
                    StructuralIdentifiability.Nemo.fmpq_mpoly,
                },
            },
        }()
    )
    eqs_expr = []
    for eq in equations
        if eq.head != :(=)
            throw("Problem with parsing at $eq")
        end
        lhs, rhs = eq.args[1:2]
        loc_all_symb = macrohelper_extract_vars([rhs])[4]
        to_insert = undef
        if lhs in x_vars
            to_insert = x_dict
        elseif lhs in y_vars
            to_insert = y_dict
        else
            throw("Unknown left-hand side $lhs")
        end

        uniqueness_check_expr = quote
            if haskey($to_insert, $lhs)
                throw(
                    DomainError(
                        $lhs,
                        "The variable occurs twice in the left-hand-side of the ODE system",
                    ),
                )
            end
        end
        push!(eqs_expr, uniqueness_check_expr)
        if isempty(loc_all_symb)
            push!(eqs_expr, :($to_insert[$lhs] = $R($rhs)))
        else
            push!(eqs_expr, :($to_insert[$lhs] = ($rhs)))
        end
    end

    for n in all_symb_no_t
        if !Base.isidentifier(n)
            throw(
                ArgumentError(
                    "The names of the variables will be injected into the global scope, so their name must be allowed Julia names, $n is not",
                ),
            )
        end
    end

    logging_exprs = [
        :(
            StructuralIdentifiability.Logging.with_logger(
                StructuralIdentifiability._si_logger[],
            ) do
                @info "Summary of the model:"
                @info "State variables: " * $(join(map(string, collect(x_vars)), ", "))
                @info "Parameters: " * $(join(map(string, collect(params)), ", "))
                @info "Inputs: " * $(join(map(string, collect(u_vars)), ", "))
                @info "Outputs: " * $(join(map(string, collect(y_vars)), ", "))
            end
        ),
    ]
    # creating the ode object
    ode_expr = :(StructuralIdentifiability.ODE{StructuralIdentifiability.Nemo.fmpq_mpoly}(
        $vx,
        $vy,
        $x_dict,
        $y_dict,
        Array{StructuralIdentifiability.Nemo.fmpq_mpoly}([$(u_vars...)]),
    ))

    result = Expr(
        :block,
        logging_exprs...,
        exp_ring,
        assignments...,
        x_var_expr,
        y_var_expr,
        x_dict_create_expr,
        y_dict_create_expr,
        eqs_expr...,
        ode_expr,
    )
    return esc(result)
end

#------------------------------------------------------------------------------

function Base.show(io::IO, ode::ODE)
    for x in ode.x_vars
        if endswith(var_to_str(x), "(t)")
            print(io, var_to_str(x)[1:(end - 3)] * "'(t) = ")
        else
            print(io, var_to_str(x) * "' = ")
        end
        print(io, ode.x_equations[x])
        print(io, "\n")
    end
    for y in ode.y_vars
        print(io, var_to_str(y) * " = ")
        print(io, ode.y_equations[y])
        print(io, "\n")
    end
end

#------------------------------------------------------------------------------
"""
    function mtk_to_si(de::ModelingToolkit.AbstractTimeDependentSystem, measured_quantities::Array{ModelingToolkit.Equation})
    function mtk_to_si(de::ModelingToolkit.AbstractTimeDependentSystem, measured_quantities::Array{SymbolicUtils.BasicSymbolic})

Input:
- `de` - ModelingToolkit.AbstractTimeDependentSystem, a system for identifiability query
- `measured_quantities` - array of output functions (as equations of just functions)

Output:
- `ODE` object containing required data for identifiability assessment
- `conversion` dictionary from the symbols in the input MTK model to the variable
  involved in the produced `ODE` object
"""
function mtk_to_si(
    de::ModelingToolkit.AbstractTimeDependentSystem,
    measured_quantities::Array{ModelingToolkit.Equation},
)
    return __mtk_to_si(
        de,
        [(replace(string(e.lhs), "(t)" => ""), e.rhs) for e in measured_quantities],
    )
end

function mtk_to_si(
    de::ModelingToolkit.AbstractTimeDependentSystem,
    measured_quantities::Array{<:Symbolics.Num},
)
    return __mtk_to_si(
        de,
        [("y$i", Symbolics.value(e)) for (i, e) in enumerate(measured_quantities)],
    )
end

function mtk_to_si(
    de::ModelingToolkit.AbstractTimeDependentSystem,
    measured_quantities::Array{<:SymbolicUtils.BasicSymbolic},
)
    return __mtk_to_si(de, [("y$i", e) for (i, e) in enumerate(measured_quantities)])
end

#------------------------------------------------------------------------------
# old name kept for compatibility purposes

function preprocess_ode(
    de::ModelingToolkit.AbstractTimeDependentSystem,
    measured_quantities::Array{ModelingToolkit.Equation},
)
    @warn "Function `preprocess_ode` has been renamed to `mtk_to_si`. The old name can be still used but will disappear in the future releases."
    return mtk_to_si(de, measured_quantities)
end

function preprocess_ode(
    de::ModelingToolkit.AbstractTimeDependentSystem,
    measured_quantities::Array{<:Symbolics.Num},
)
    @warn "Function `preprocess_ode` has been renamed to `mtk_to_si`. The old name can be still used but will disappear in the future releases."
    return mtk_to_si(de, measured_quantities)
end

function preprocess_ode(
    de::ModelingToolkit.AbstractTimeDependentSystem,
    measured_quantities::Array{<:SymbolicUtils.BasicSymbolic},
)
    @warn "Function `preprocess_ode` has been renamed to `mtk_to_si`. The old name can be still used but will disappear in the future releases."
    return mtk_to_si(de, measured_quantities)
end

#------------------------------------------------------------------------------
"""
    function __mtk_to_si(de::ModelingToolkit.AbstractTimeDependentSystem, measured_quantities::Array{Tuple{String, SymbolicUtils.BasicSymbolic}})

Input:
- `de` - ModelingToolkit.AbstractTimeDependentSystem, a system for identifiability query
- `measured_quantities` - array of input function in the form (name, expression)

Output:
- `ODE` object containing required data for identifiability assessment
- `conversion` dictionary from the symbols in the input MTK model to the variable
  involved in the produced `ODE` object
"""
function __mtk_to_si(
    de::ModelingToolkit.AbstractTimeDependentSystem,
    measured_quantities::Array{<:Tuple{String, <:SymbolicUtils.BasicSymbolic}},
)
    polytype = StructuralIdentifiability.Nemo.fmpq_mpoly
    fractype = StructuralIdentifiability.Nemo.Generic.Frac{polytype}
    diff_eqs =
        filter(eq -> !(ModelingToolkit.isoutput(eq.lhs)), ModelingToolkit.equations(de))
    # performing full structural simplification
    if length(observed(de)) > 0
        rules = Dict(s.lhs => s.rhs for s in observed(de))
        while any([
            length(intersect(get_variables(r), keys(rules))) > 0 for r in values(rules)
        ])
            rules = Dict(k => SymbolicUtils.substitute(v, rules) for (k, v) in rules)
        end
        diff_eqs = [SymbolicUtils.substitute(eq, rules) for eq in diff_eqs]
    end

    y_functions = [each[2] for each in measured_quantities]
    inputs = filter(v -> ModelingToolkit.isinput(v), ModelingToolkit.states(de))
    state_vars = filter(
        s -> !(ModelingToolkit.isinput(s) || ModelingToolkit.isoutput(s)),
        ModelingToolkit.states(de),
    )
    params = ModelingToolkit.parameters(de)
    t = ModelingToolkit.arguments(diff_eqs[1].lhs)[1]
    params_from_measured_quantities = union(
        [filter(s -> !istree(s), get_variables(y[2])) for y in measured_quantities]...,
    )
    params = union(params, params_from_measured_quantities)

    input_symbols = vcat(state_vars, inputs, params)
    generators = vcat(string.(input_symbols), [e[1] for e in measured_quantities])
    generators = map(g -> replace(g, "(t)" => ""), generators)
    R, gens_ = Nemo.PolynomialRing(Nemo.QQ, generators)
    y_vars = Vector{polytype}([str_to_var(e[1], R) for e in measured_quantities])
    symb2gens = Dict(input_symbols .=> gens_[1:length(input_symbols)])

    x_vars = Vector{polytype}()

    state_eqn_dict = Dict{polytype, Union{polytype, fractype}}()
    out_eqn_dict = Dict{polytype, Union{polytype, fractype}}()

    for i in 1:length(diff_eqs)
        x = substitute(state_vars[i], symb2gens)
        push!(x_vars, x)
        if !(diff_eqs[i].rhs isa Number)
            state_eqn_dict[x] = eval_at_nemo(diff_eqs[i].rhs, symb2gens)
        else
            state_eqn_dict[x] = R(diff_eqs[i].rhs)
        end
    end
    for i in 1:length(measured_quantities)
        out_eqn_dict[y_vars[i]] = eval_at_nemo(measured_quantities[i][2], symb2gens)
    end

    inputs_ = [substitute(each, symb2gens) for each in inputs]
    if isequal(length(inputs_), 0)
        inputs_ = Vector{polytype}()
    end
    return (
        StructuralIdentifiability.ODE{polytype}(
            x_vars,
            y_vars,
            state_eqn_dict,
            out_eqn_dict,
            inputs_,
        ),
        symb2gens,
    )
end
