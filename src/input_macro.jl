function _extract_aux!(funcs, all_symb, eq; ders_ok = false, type = :ode)
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
                if type != :ode
                    throw(
                        Base.ArgumentError(
                            "Derivative are not expected in the discrete case",
                        ),
                    )
                end
                push!(all_symb, f)
            elseif @capture(x, f_(t + 1))
                if !ders_ok
                    throw(Base.ArgumentError("Shifts are not allowed in the right-hand side"))
                end
                if type != :dds
                    throw(
                        Base.ArgumentError(
                            "Shifts are not expected in the differential case",
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
  For an expression of the form f'(t)/f(t + 1) or f(t) returns (f, true) and (f, false), resp
"""
function _get_var(expr, type = :ode)
    if @capture(expr, f_'(t))
        @assert type == :ode
        return (f, true)
    end
    if @capture(expr, f_(t + 1))
        @assert type == :dds
        return (f, true)
    end
    if @capture(expr, f_(t))
        return (f, false)
    end
    error("cannot extract the single function name from $expr")
end

function macrohelper_extract_vars(equations::Array{Expr, 1}, type = :ode)
    funcs, all_symb = Set(), Set()
    x_vars, y_vars = Vector(), Vector()
    for eq in equations
        if eq.head != :(=)
            _extract_aux!(funcs, all_symb, eq, type = type)
        else
            lhs, rhs = eq.args[1:2]
            _extract_aux!(funcs, all_symb, lhs, ders_ok = true, type = type)
            _extract_aux!(funcs, all_symb, rhs, type = type)
            (v, is_state) = _get_var(lhs, type)
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

function macrohelper_extract_vars(equations::Array{Symbol, 1}, type = :ode)
    return macrohelper_extract_vars(map(Expr, equations), type)
end

#------------------------------------------------------------------------------

function macrohelper_clean(ex::Expr)
    ex = MacroTools.postwalk(x -> @capture(x, f_'(t)) ? f : x, ex)
    ex = MacroTools.postwalk(x -> @capture(x, f_(t + 1)) ? f : x, ex)
    ex = MacroTools.postwalk(x -> @capture(x, f_(t)) ? f : x, ex)
    ex = MacroTools.postwalk(x -> x == :(/) ? :(//) : x, ex)
    ex = MacroTools.postwalk(x -> x isa Float64 ? rationalize(x) : x, ex)
    return ex
end

#------------------------------------------------------------------------------

function generate_model_code(type, ex::Expr...)
    @assert type in (:ode, :dds)
    equations = [ex...]
    x_vars, y_vars, u_vars, all_symb = macrohelper_extract_vars(equations, type)
    time_dependent = vcat(x_vars, y_vars, u_vars)
    params = sort([s for s in all_symb if !(s in time_dependent)])
    all_symb_no_t = vcat(time_dependent, params)
    all_symb_with_t = vcat([:($s(t)) for s in time_dependent], params)

    # creating the polynomial ring
    R = gensym()
    vars_aux = gensym()
    exp_ring = :(
        ($R, $vars_aux) = StructuralIdentifiability.Nemo.polynomial_ring(
            StructuralIdentifiability.Nemo.QQ,
            map(string, $all_symb_with_t),
        )
    )
    assignments = [:($(all_symb_no_t[i]) = $vars_aux[$i]) for i in 1:length(all_symb_no_t)]

    # setting x_vars and y_vars in the right order
    vx = gensym()
    vy = gensym()
    x_var_expr =
        :($vx = Vector{StructuralIdentifiability.Nemo.QQMPolyRingElem}([$(x_vars...)]))
    y_var_expr =
        :($vy = Vector{StructuralIdentifiability.Nemo.QQMPolyRingElem}([$(y_vars...)]))

    # preparing equations
    equations = map(macrohelper_clean, equations)
    x_dict = gensym()
    y_dict = gensym()
    x_dict_create_expr = :(
        $x_dict = Dict{
            StructuralIdentifiability.Nemo.QQMPolyRingElem,
            Union{
                StructuralIdentifiability.Nemo.QQMPolyRingElem,
                StructuralIdentifiability.AbstractAlgebra.Generic.Frac{
                    StructuralIdentifiability.Nemo.QQMPolyRingElem,
                },
            },
        }()
    )
    y_dict_create_expr = :(
        $y_dict = Dict{
            StructuralIdentifiability.Nemo.QQMPolyRingElem,
            Union{
                StructuralIdentifiability.Nemo.QQMPolyRingElem,
                StructuralIdentifiability.AbstractAlgebra.Generic.Frac{
                    StructuralIdentifiability.Nemo.QQMPolyRingElem,
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
        loc_all_symb = macrohelper_extract_vars([rhs], type)[4]
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
    # creating the ode/dds object
    obj_type = Dict(:ode => :ODE, :dds => :DDS)
    ds_expr = :(StructuralIdentifiability.$(obj_type[type]){
        StructuralIdentifiability.Nemo.QQMPolyRingElem,
    }(
        $vx,
        $vy,
        $x_dict,
        $y_dict,
        Array{StructuralIdentifiability.Nemo.QQMPolyRingElem}([$(u_vars...)]),
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
        ds_expr,
    )
    return result
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
- `a`, `b`, `c` are time-independent parameters

"""
macro ODEmodel(ex::Expr...)
    return esc(generate_model_code(:ode, ex...))
end

#------------------------------------------------------------------------------

"""
    macro DDSmodel

Macro for creating a DDS (discrete dynamical system) 
from a list of equations.
It also injects all variables into the global scope.

## Example

Creating a simple `DDS`:

```jldoctest
using StructuralIdentifiability

dds = @DDSmodel(
    x1(t + 1) = a * x1(t) + u(t),
    x2(t + 1) = b * x2(t) + c*x1(t)*x2(t),
    y(t) = x1(t)
)
```

Here,
- `x1`, `x2` are state variables
- `y` is an output variable
- `u` is an input variable
- `a`, `b`, `c` are time-indepdendent parameters

"""
macro DDSmodel(ex::Expr...)
    return esc(generate_model_code(:dds, ex...))
end

#------------------------------------------------------------------------------
