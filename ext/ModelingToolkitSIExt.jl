module ModelingToolkitSIExt

using DataStructures
using Logging
using Nemo
using Random
using StructuralIdentifiability
using StructuralIdentifiability: str_to_var, parent_ring_change, eval_at_dict
using StructuralIdentifiability: restart_logging, _si_logger, reset_timings, _to
using TimerOutputs

if isdefined(Base, :get_extension)
    using ModelingToolkit
else
    using ..ModelingToolkit
end

# ------------------------------------------------------------------------------

function StructuralIdentifiability.eval_at_nemo(e::Num, vals::Dict)
    e = Symbolics.value(e)
    return eval_at_nemo(e, vals)
end

function StructuralIdentifiability.eval_at_nemo(e::SymbolicUtils.BasicSymbolic, vals::Dict)
    if Symbolics.iscall(e)
        # checking if it is a function of the form x(t), a bit dirty
        if length(Symbolics.arguments(e)) == 1 && "$(first(Symbolics.arguments(e)))" == "t"
            return vals[e]
        end
        # checking if this is a vector entry like x(t)[1]
        if Symbolics.operation(e) == getindex
            return vals[e]
        end
        # otherwise, this is a term
        args = map(
            a -> StructuralIdentifiability.eval_at_nemo(a, vals),
            Symbolics.arguments(e),
        )
        if Symbolics.operation(e) in (+, -, *)
            return Symbolics.operation(e)(args...)
        elseif isequal(Symbolics.operation(e), /)
            return //(args...)
        elseif isequal(Symbolics.operation(e), ^)
            if args[2] >= 0
                return args[1]^args[2]
            end
            return 1 // args[1]^(-args[2])
        end
        throw(Base.ArgumentError("Function $(Symbolics.operation(e)) is not supported"))
    elseif e isa Symbolics.Symbolic
        return get(vals, e, e)
    end
end

function StructuralIdentifiability.eval_at_nemo(e::Union{Integer, Rational}, vals::Dict)
    return e
end

function StructuralIdentifiability.eval_at_nemo(
    e::Union{Float16, Float32, Float64},
    vals::Dict,
)
    if isequal(e % 1, 0)
        out = Int(e)
    else
        out = rationalize(e)
    end
    @warn "Floating point value $e will be converted to $(out)."
    return out
end

function get_measured_quantities(ode::ModelingToolkit.ODESystem)
    if any(ModelingToolkit.isoutput(eq.lhs) for eq in ModelingToolkit.equations(ode))
        @info "Measured quantities are not provided, trying to find the outputs in input ODE."
        return filter(
            eq -> (ModelingToolkit.isoutput(eq.lhs)),
            ModelingToolkit.equations(ode),
        )
    else
        throw(
            error(
                "Measured quantities (output functions) were not provided and no outputs were found.",
            ),
        )
    end
end

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
function StructuralIdentifiability.mtk_to_si(
    de::ModelingToolkit.AbstractTimeDependentSystem,
    measured_quantities::Array{ModelingToolkit.Equation},
)
    return __mtk_to_si(
        de,
        [(replace(string(e.lhs), "(t)" => ""), e.rhs) for e in measured_quantities],
    )
end

function StructuralIdentifiability.mtk_to_si(
    de::ModelingToolkit.AbstractTimeDependentSystem,
    measured_quantities::Array{<:Symbolics.Num},
)
    return __mtk_to_si(
        de,
        [("y$i", Symbolics.value(e)) for (i, e) in enumerate(measured_quantities)],
    )
end

function StructuralIdentifiability.mtk_to_si(
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
- `measured_quantities` - array of output function in the form (name, expression)

Output:
- `ODE` object containing required data for identifiability assessment
- `conversion` dictionary from the symbols in the input MTK model to the variable
  involved in the produced `ODE` object
"""
function __mtk_to_si(
    de::ModelingToolkit.AbstractTimeDependentSystem,
    measured_quantities::Array{<:Tuple{String, <:SymbolicUtils.BasicSymbolic}},
)
    polytype = StructuralIdentifiability.Nemo.QQMPolyRingElem
    fractype = StructuralIdentifiability.Nemo.Generic.FracFieldElem{polytype}
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
    inputs = filter(v -> ModelingToolkit.isinput(v), ModelingToolkit.unknowns(de))
    state_vars = filter(
        s -> !(ModelingToolkit.isinput(s) || ModelingToolkit.isoutput(s)),
        ModelingToolkit.unknowns(de),
    )
    params = ModelingToolkit.parameters(de)
    t = ModelingToolkit.arguments(diff_eqs[1].lhs)[1]
    params_from_measured_quantities = union(
        [filter(s -> !iscall(s), get_variables(y[2])) for y in measured_quantities]...,
    )
    params = union(params, params_from_measured_quantities)

    input_symbols = vcat(state_vars, inputs, params)
    generators = vcat(string.(input_symbols), [e[1] for e in measured_quantities])
    generators = map(g -> replace(g, "(t)" => ""), generators)
    R, gens_ = Nemo.polynomial_ring(Nemo.QQ, generators)
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
# -----------------------------------------------------------------------------
"""
    function assess_local_identifiability(ode::ModelingToolkit.ODESystem; measured_quantities=Array{ModelingToolkit.Equation}[], funcs_to_check=Array{}[], prob_threshold::Float64=0.99, type=:SE, loglevel=Logging.Info)

Input:
- `ode` - the ODESystem object from ModelingToolkit
- `measured_quantities` - the measurable outputs of the model
- `funcs_to_check` - functions of parameters for which to check identifiability
- `prob_threshold` - probability of correctness
- `type` - identifiability type (`:SE` for single-experiment, `:ME` for multi-experiment)
- `loglevel` - the minimal level of log messages to display (`Logging.Info` by default)

Output:
- for `type=:SE`, the result is an (ordered) dictionary from each parameter to boolean;
- for `type=:ME`, the result is a tuple with the dictionary as in `:SE` case and array of number of experiments.

The function determines local identifiability of parameters in `funcs_to_check` or all possible parameters if `funcs_to_check` is empty

The result is correct with probability at least `prob_threshold`.

`type` can be either `:SE` (single-experiment identifiability) or `:ME` (multi-experiment identifiability).
The return value is a tuple consisting of the array of bools and the number of experiments to be performed.
"""
function StructuralIdentifiability.assess_local_identifiability(
    ode::ModelingToolkit.ODESystem;
    measured_quantities = Array{ModelingToolkit.Equation}[],
    funcs_to_check = Array{}[],
    prob_threshold::Float64 = 0.99,
    type = :SE,
    loglevel = Logging.Info,
)
    restart_logging(loglevel = loglevel)
    with_logger(_si_logger[]) do
        return _assess_local_identifiability(
            ode,
            measured_quantities = measured_quantities,
            funcs_to_check = funcs_to_check,
            prob_threshold = prob_threshold,
            type = type,
        )
    end
end

@timeit _to function _assess_local_identifiability(
    ode::ModelingToolkit.ODESystem;
    measured_quantities = Array{ModelingToolkit.Equation}[],
    funcs_to_check = Array{}[],
    prob_threshold::Float64 = 0.99,
    type = :SE,
)
    if length(measured_quantities) == 0
        if any(ModelingToolkit.isoutput(eq.lhs) for eq in ModelingToolkit.equations(ode))
            @info "Measured quantities are not provided, trying to find the outputs in input ODE."
            measured_quantities = filter(
                eq -> (ModelingToolkit.isoutput(eq.lhs)),
                ModelingToolkit.equations(ode),
            )
        else
            throw(
                error(
                    "Measured quantities (output functions) were not provided and no outputs were found.",
                ),
            )
        end
    end
    if length(funcs_to_check) == 0
        funcs_to_check = vcat(
            [e for e in ModelingToolkit.unknowns(ode) if !ModelingToolkit.isoutput(e)],
            ModelingToolkit.parameters(ode),
        )
    end
    ode, conversion = mtk_to_si(ode, measured_quantities)
    funcs_to_check_ = [eval_at_nemo(x, conversion) for x in funcs_to_check]

    if isequal(type, :SE)
        result = StructuralIdentifiability._assess_local_identifiability(
            ode,
            funcs_to_check = funcs_to_check_,
            prob_threshold = prob_threshold,
            type = type,
        )
        nemo2mtk = Dict(funcs_to_check_ .=> funcs_to_check)
        out_dict =
            OrderedDict(nemo2mtk[param] => result[param] for param in funcs_to_check_)
        return out_dict
    elseif isequal(type, :ME)
        result, bd = StructuralIdentifiability._assess_local_identifiability(
            ode,
            funcs_to_check = funcs_to_check_,
            prob_threshold = prob_threshold,
            type = type,
        )
        nemo2mtk = Dict(funcs_to_check_ .=> funcs_to_check)
        out_dict =
            OrderedDict(nemo2mtk[param] => result[param] for param in funcs_to_check_)
        return (out_dict, bd)
    end
end

# ------------------------------------------------------------------------------

"""
    assess_identifiability(ode::ModelingToolkit.ODESystem; measured_quantities=Array{ModelingToolkit.Equation}[], funcs_to_check=[], known_ic=[], prob_threshold = 0.99, loglevel=Logging.Info)

Input:
- `ode` - the ModelingToolkit.ODESystem object that defines the model
- `measured_quantities` - the output functions of the model
- `funcs_to_check` - functions of parameters for which to check the identifiability
- `known_ic` - functions, for which initial conditions are assumed to be known
- `prob_threshold` - probability of correctness.
- `loglevel` - the minimal level of log messages to display (`Logging.Info` by default)

Assesses identifiability (both local and global) of a given ODE model (parameters detected automatically). The result is guaranteed to be correct with the probability
at least `prob_threshold`.
If known initial conditions are provided, the identifiability results for the states will also hold at `t = 0`
"""
function StructuralIdentifiability.assess_identifiability(
    ode::ModelingToolkit.ODESystem;
    measured_quantities = Array{ModelingToolkit.Equation}[],
    funcs_to_check = [],
    known_ic = [],
    prob_threshold = 0.99,
    loglevel = Logging.Info,
)
    restart_logging(loglevel = loglevel)
    with_logger(_si_logger[]) do
        return _assess_identifiability(
            ode,
            measured_quantities = measured_quantities,
            funcs_to_check = funcs_to_check,
            known_ic = known_ic,
            prob_threshold = prob_threshold,
        )
    end
end

function _assess_identifiability(
    ode::ModelingToolkit.ODESystem;
    measured_quantities = Array{ModelingToolkit.Equation}[],
    funcs_to_check = [],
    known_ic = [],
    prob_threshold = 0.99,
)
    if isempty(measured_quantities)
        measured_quantities = get_measured_quantities(ode)
    end

    ode, conversion = mtk_to_si(ode, measured_quantities)
    conversion_back = Dict(v => k for (k, v) in conversion)
    if isempty(funcs_to_check)
        funcs_to_check = [conversion_back[x] for x in [ode.x_vars..., ode.parameters...]]
    end
    funcs_to_check_ = [eval_at_nemo(each, conversion) for each in funcs_to_check]

    known_ic_ = [eval_at_nemo(each, conversion) for each in known_ic]

    nemo2mtk = Dict(funcs_to_check_ .=> funcs_to_check)
    result = nothing
    if isempty(known_ic)
        result = StructuralIdentifiability._assess_identifiability(
            ode,
            funcs_to_check = funcs_to_check_,
            prob_threshold = prob_threshold,
        )
        return OrderedDict(nemo2mtk[param] => result[param] for param in funcs_to_check_)
    else
        result = StructuralIdentifiability._assess_identifiability_kic(
            ode,
            known_ic_,
            funcs_to_check = funcs_to_check_,
            prob_threshold = prob_threshold,
        )
    end
    nemo2mtk = Dict(funcs_to_check_ .=> funcs_to_check)
    out_dict = OrderedDict(nemo2mtk[param] => result[param] for param in funcs_to_check_)
    if length(known_ic) > 0
        @warn "Since known initial conditions were provided, identifiability of states (e.g., `x(t)`) is at t = 0 only !"
        t = SymbolicUtils.Sym{Real}(:t)
        out_dict = OrderedDict(substitute(k, Dict(t => 0)) => v for (k, v) in out_dict)
    end
    return out_dict
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
- `dds` - the DiscreteSystem object from ModelingToolkit
- `measured_quantities` - the measurable outputs of the model
- `funcs_to_check` - functions of parameters for which to check identifiability (all parameters and states if not specified)
- `known_ic` - functions (of states and parameter) whose initial conditions are assumed to be known
- `prob_threshold` - probability of correctness

Output:
- the result is an (ordered) dictionary from each function to to boolean;

The result is correct with probability at least `prob_threshold`.
"""
function StructuralIdentifiability.assess_local_identifiability(
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

    dds_aux_ode, conversion = mtk_to_si(dds, measured_quantities)
    dds_aux = StructuralIdentifiability.DDS{QQMPolyRingElem}(dds_aux_ode)
    if length(funcs_to_check) == 0
        params = parameters(dds)
        params_from_measured_quantities = union(
            [filter(s -> !iscall(s), get_variables(y)) for y in measured_quantities]...,
        )
        funcs_to_check = vcat(
            [
                x for x in unknowns(dds) if
                conversion[x] in StructuralIdentifiability.x_vars(dds_aux)
            ],
            union(params, params_from_measured_quantities),
        )
    end
    funcs_to_check_ = [eval_at_nemo(x, conversion) for x in funcs_to_check]
    known_ic_ = [eval_at_nemo(x, conversion) for x in known_ic]

    result = StructuralIdentifiability._assess_local_identifiability_discrete_aux(
        dds_aux,
        funcs_to_check_,
        known_ic_,
        prob_threshold,
    )
    nemo2mtk = Dict(funcs_to_check_ .=> funcs_to_check)
    out_dict = OrderedDict(nemo2mtk[param] => result[param] for param in funcs_to_check_)
    if length(known_ic) > 0
        @warn "Since known initial conditions were provided, identifiability of states (e.g., `x(t)`) is at t = 0 only !"
        t = SymbolicUtils.Sym{Real}(:t)
        out_dict = OrderedDict(substitute(k, Dict(t => 0)) => v for (k, v) in out_dict)
    end
    return out_dict
end

# ------------------------------------------------------------------------------

"""
    find_identifiable_functions(ode::ModelingToolkit.ODESystem; measured_quantities=[], known_ic=[], options...)

Finds all functions of parameters/states that are identifiable in the given ODE
system.

## Options

This functions takes the following optional arguments:
- `measured_quantities` - the output functions of the model.
- `known_ic` - a list of functions whose initial conditions are assumed to be known,
  then the returned identifiable functions will be functions of parameters and
  initial conditions, not states (this is an experimental functionality).
- `loglevel` - the verbosity of the logging
  (can be Logging.Error, Logging.Warn, Logging.Info, Logging.Debug)

## Example

```jldoctest
using StructuralIdentifiability
using ModelingToolkit

@parameters a01 a21 a12
@variables t x0(t) x1(t) y1(t) [output = true]
D = Differential(t)

eqs = [
    D(x0) ~ -(a01 + a21) * x0 + a12 * x1,
    D(x1) ~ a21 * x0 - a12 * x1, y1 ~ x0
]
de = ODESystem(eqs, t, name = :Test)

find_identifiable_functions(de, measured_quantities = [y1 ~ x0])

# prints
2-element Vector{Num}:
         a01*a12
 a01 + a12 + a21
```
"""
function StructuralIdentifiability.find_identifiable_functions(
    ode::ModelingToolkit.ODESystem;
    measured_quantities = Array{ModelingToolkit.Equation}[],
    known_ic = [],
    prob_threshold::Float64 = 0.99,
    seed = 42,
    with_states = false,
    simplify = :standard,
    rational_interpolator = :VanDerHoevenLecerf,
    loglevel = Logging.Info,
)
    restart_logging(loglevel = loglevel)
    reset_timings()
    with_logger(_si_logger[]) do
        return _find_identifiable_functions(
            ode,
            measured_quantities = measured_quantities,
            known_ic = known_ic,
            prob_threshold = prob_threshold,
            seed = seed,
            with_states = with_states,
            simplify = simplify,
            rational_interpolator = rational_interpolator,
        )
    end
end

function _find_identifiable_functions(
    ode::ModelingToolkit.ODESystem;
    measured_quantities = Array{ModelingToolkit.Equation}[],
    known_ic = Array{Symbolics.Num}[],
    prob_threshold::Float64 = 0.99,
    seed = 42,
    with_states = false,
    simplify = :standard,
    rational_interpolator = :VanDerHoevenLecerf,
)
    Random.seed!(seed)
    if isempty(measured_quantities)
        measured_quantities = get_measured_quantities(ode)
    end
    ode, conversion = mtk_to_si(ode, measured_quantities)
    known_ic_ = [eval_at_nemo(each, conversion) for each in known_ic]
    result = nothing
    if isempty(known_ic)
        result = StructuralIdentifiability._find_identifiable_functions(
            ode,
            simplify = simplify,
            prob_threshold = prob_threshold,
            seed = seed,
            with_states = with_states,
            rational_interpolator = rational_interpolator,
        )
    else
        result = StructuralIdentifiability._find_identifiable_functions_kic(
            ode,
            known_ic_,
            simplify = simplify,
            prob_threshold = prob_threshold,
            seed = seed,
            rational_interpolator = rational_interpolator,
        )
    end
    result = [parent_ring_change(f, ode.poly_ring) for f in result]
    nemo2mtk = Dict(v => Num(k) for (k, v) in conversion)
    out_funcs = [eval_at_dict(func, nemo2mtk) for func in result]
    if length(known_ic) > 0
        @warn "Since known initial conditions were provided, identifiability of states (e.g., `x(t)`) is at t = 0 only !"
        t = SymbolicUtils.Sym{Real}(:t)
        out_funcs = [substitute(f, Dict(t => 0)) for f in out_funcs]
    end

    return out_funcs
end

end
