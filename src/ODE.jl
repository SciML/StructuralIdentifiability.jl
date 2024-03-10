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
    x_equations::Dict{P, <:ExtendedFraction{P}}
    y_equations::Dict{P, <:ExtendedFraction{P}}

    function ODE{P}(
        x_vars::Array{P, 1},
        y_vars::Array{P, 1},
        x_eqs::Dict{P, <:ExtendedFraction{P}},
        y_eqs::Dict{P, <:ExtendedFraction{P}},
        inputs::Array{P, 1},
    ) where {P <: MPolyRingElem{<:FieldElem}}
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
        x_eqs::Dict{P, <:ExtendedFraction{P}},
        y_eqs::Dict{P, <:ExtendedFraction{P}},
        inputs::Array{P, 1},
    ) where {P <: MPolyRingElem{<:FieldElem}}
        x_vars = collect(keys(x_eqs))
        y_vars = collect(keys(y_eqs))
        return ODE{P}(x_vars, y_vars, x_eqs, y_eqs, inputs)
    end
end

function Base.parent(ode::ODE)
    return ode.poly_ring
end

#------------------------------------------------------------------------------

function add_outputs(
    ode::ODE{P},
    extra_y::Dict{String, <:RingElem},
) where {P <: MPolyRingElem}
    new_var_names =
        vcat(collect(map(var_to_str, gens(ode.poly_ring))), collect(keys(extra_y)))
    new_ring, new_vars = Nemo.polynomial_ring(base_ring(ode.poly_ring), new_var_names)

    new_x = Array{P, 1}([parent_ring_change(x, new_ring) for x in ode.x_vars])
    new_x_eqs = Dict{P, ExtendedFraction{P}}(
        parent_ring_change(x, new_ring) => parent_ring_change(f, new_ring) for
        (x, f) in ode.x_equations
    )
    new_y = Array{P, 1}([parent_ring_change(y, new_ring) for y in ode.y_vars])
    for y in keys(extra_y)
        push!(new_y, str_to_var(y, new_ring))
    end
    new_y_eqs = Dict{P, ExtendedFraction{P}}(
        parent_ring_change(y, new_ring) => parent_ring_change(g, new_ring) for
        (y, g) in ode.y_equations
    )
    extra_y_eqs = Dict{P, ExtendedFraction{P}}(
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
) where {T <: FieldElem, P <: MPolyRingElem{T}}
    new_vars =
        map(var_to_str, [v for v in gens(ode.poly_ring) if !(v in keys(param_values))])
    small_ring, small_vars = Nemo.polynomial_ring(base_ring(ode.poly_ring), new_vars)
    eval_dict =
        Dict(str_to_var(v, ode.poly_ring) => str_to_var(v, small_ring) for v in new_vars)
    merge!(eval_dict, Dict(p => small_ring(val) for (p, val) in param_values))

    return ODE{P}(
        Dict{P, ExtendedFraction{P}}(
            eval_at_dict(v, eval_dict) => eval_at_dict(f, eval_dict) for
            (v, f) in ode.x_equations
        ),
        Dict{P, ExtendedFraction{P}}(
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
) where {P <: MPolyRingElem}
    new_values = Dict{P, QQFieldElem}(x => _to_rational(v) for (x, v) in param_values)
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
) where {T <: FieldElem, P <: MPolyRingElem{T}}
    new_varnames = map(var_to_str, vcat(ode.x_vars, ode.u_vars))
    append!(new_varnames, map(v -> var_to_str(v) * "_dot", ode.x_vars))

    new_ring, new_vars = Nemo.polynomial_ring(base_ring(ode.poly_ring), new_varnames)
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
) where {P <: MPolyRingElem{<:FieldElem}}
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
function _reduce_mod_p(poly::QQMPolyRingElem, p::Int)
    den = denominator(poly)
    num = change_base_ring(Nemo.ZZ, den * poly)
    if Nemo.Native.GF(p)(den) == 0
        throw(Base.ArgumentError("Prime $p divides the denominator of $poly"))
    end
    return change_base_ring(Nemo.Native.GF(p), num) * (1 // Nemo.Native.GF(p)(den))
end

function _reduce_mod_p(rat::Generic.FracFieldElem{QQMPolyRingElem}, p::Int)
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
function reduce_ode_mod_p(ode::ODE{<:MPolyRingElem{Nemo.QQFieldElem}}, p::Int)
    new_ring, new_vars =
        Nemo.polynomial_ring(Nemo.Native.GF(p), map(var_to_str, gens(ode.poly_ring)))
    new_type = typeof(new_vars[1])
    new_inputs = map(u -> switch_ring(u, new_ring), ode.u_vars)
    new_x = map(x -> switch_ring(x, new_ring), ode.x_vars)
    new_y = map(y -> switch_ring(y, new_ring), ode.y_vars)
    new_x_eqs = Dict{new_type, ExtendedFraction{new_type}}()
    new_y_eqs = Dict{new_type, ExtendedFraction{new_type}}()
    for (old, new) in Dict(ode.x_equations => new_x_eqs, ode.y_equations => new_y_eqs)
        for (v, f) in old
            new_v = switch_ring(v, new_ring)
            new[new_v] = _reduce_mod_p(f, p)
        end
    end
    return ODE{new_type}(new_x, new_y, new_x_eqs, new_y_eqs, new_inputs)
end

#------------------------------------------------------------------------------

function Base.show(io::IO, ode::ODE)
    for x in ode.x_vars
        if endswith(var_to_str(x), "(t)")
            print(io, chopsuffix(var_to_str(x), "(t)") * "'(t) = ")
            # print(io, var_to_str(x)[1:(end - 3)] * "'(t) = ")
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
