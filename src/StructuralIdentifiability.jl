module StructuralIdentifiability

# General purpose packages
using Base.Iterators
using Combinatorics
using DataStructures
using IterTools
using LinearAlgebra
using Logging
using MacroTools
using Primes
using Random
using TimerOutputs

# Algebra packages
using AbstractAlgebra
using Nemo
using Groebner
using ParamPunPam
using ParamPunPam: reduce_mod_p!, specialize_mod_p, AbstractBlackboxIdeal
ParamPunPam.enable_progressbar(false)

using ModelingToolkit

# defining a model
export ODE, @ODEmodel, mtk_to_si

# assessing identifiability
export assess_local_identifiability, assess_identifiability

# auxuliary function
export set_parameter_values

# extra functionality
export find_ioequations, find_identifiable_functions

# exporting to other formats
export print_for_maple, print_for_DAISY, print_for_COMBOS, print_for_GenSSI

# function for creating linear compartment models
export linear_compartment_model

# structure for storing the io-equations for doing reductions further
export PBRepresentation, diffreduce, io_switch!, pseudodivision

# functions for finding submodels of a given model and visualization of the submodels
export find_submodels

# finding identifiabile reparametrizations
export reparametrize_global

include("logging.jl")
include("util.jl")
include("power_series_utils.jl")
include("ODE.jl")
include("ODEexport.jl")
include("local_identifiability.jl")
include("wronskian.jl")
include("elimination.jl")
include("primality_check.jl")
include("io_equation.jl")
include("states.jl")
include("RationalFunctionFields/util.jl")
include("RationalFunctionFields/IdealMQS.jl")
include("RationalFunctionFields/rankings.jl")
include("RationalFunctionFields/normalforms.jl")
include("RationalFunctionFields/RationalFunctionField.jl")
include("global_identifiability.jl")
include("identifiable_functions.jl")
include("parametrizations.jl")
include("lincomp.jl")
include("pb_representation.jl")
include("submodels.jl")
include("discrete.jl")

function __init__()
    _si_logger[] = @static if VERSION >= v"1.7.0"
        Logging.ConsoleLogger(Logging.Info, show_limited = false)
    else
        Logging.ConsoleLogger(stderr, Logging.Info)
    end
end

"""
    assess_identifiability(ode; funcs_to_check = [], prob_threshold=0.99, loglevel=Logging.Info)

Input:
- `ode` - the ODE model
- `funcs_to_check` - list of functions to check identifiability for; if empty, all parameters
   and states are taken
- `prob_threshold` - probability of correctness.
- `loglevel` - the minimal level of log messages to display (`Logging.Info` by default)

Assesses identifiability of a given ODE model. The result is guaranteed to be correct with the probability
at least `prob_threshold`.
The function returns an (ordered) dictionary from the functions to check to their identifiability properties 
(one of `:nonidentifiable`, `:locally`, `:globally`).
"""
function assess_identifiability(
    ode::ODE{P};
    funcs_to_check = Vector(),
    prob_threshold::Float64 = 0.99,
    loglevel = Logging.Info,
) where {P <: MPolyElem{fmpq}}
    restart_logging(loglevel = loglevel)
    reset_timings()
    with_logger(_si_logger[]) do
        return _assess_identifiability(
            ode,
            funcs_to_check = funcs_to_check,
            prob_threshold = prob_threshold,
        )
    end
end

function _assess_identifiability(
    ode::ODE{P};
    funcs_to_check = Vector(),
    prob_threshold::Float64 = 0.99,
) where {P <: MPolyElem{fmpq}}
    p_glob = 1 - (1 - prob_threshold) * 0.9
    p_loc = 1 - (1 - prob_threshold) * 0.1

    if isempty(funcs_to_check)
        funcs_to_check = vcat(ode.x_vars, ode.parameters)
    end

    @info "Assessing local identifiability"
    trbasis = Array{fmpq_mpoly, 1}()
    runtime = @elapsed local_result = _assess_local_identifiability(
        ode,
        funcs_to_check = funcs_to_check,
        prob_threshold = p_loc,
        type = :SE,
        trbasis = trbasis,
    )
    @debug "Local identifiability assessed in $runtime seconds"
    @debug "Trasncendence basis to be specialized is $trbasis"
    # _runtime_logger[:loc_time] = runtime

    loc_id = [local_result[each] for each in funcs_to_check]
    locally_identifiable = Array{Any, 1}()
    for (loc, f) in zip(loc_id, funcs_to_check)
        if loc
            push!(locally_identifiable, f)
        end
    end

    @info "Assessing global identifiability"
    runtime = @elapsed global_result =
        assess_global_identifiability(ode, locally_identifiable, trbasis, p_glob)
    @info "Global identifiability assessed in $runtime seconds"
    _runtime_logger[:glob_time] = runtime

    result = OrderedDict{Any, Symbol}()
    glob_ind = 1
    for i in 1:length(funcs_to_check)
        if !local_result[funcs_to_check[i]]
            result[funcs_to_check[i]] = :nonidentifiable
        else
            if global_result[glob_ind]
                result[funcs_to_check[i]] = :globally
            else
                result[funcs_to_check[i]] = :locally
            end
            glob_ind += 1
        end
    end

    return result
end

"""
    assess_identifiability(ode::ModelingToolkit.ODESystem; measured_quantities=Array{ModelingToolkit.Equation}[], funcs_to_check=[], prob_threshold = 0.99, loglevel=Logging.Info)

Input:
- `ode` - the ModelingToolkit.ODESystem object that defines the model
- `measured_quantities` - the output functions of the model
- `funcs_to_check` - functions of parameters for which to check the identifiability
- `prob_threshold` - probability of correctness.
- `loglevel` - the minimal level of log messages to display (`Logging.Info` by default)

Assesses identifiability (both local and global) of a given ODE model (parameters detected automatically). The result is guaranteed to be correct with the probability
at least `prob_threshold`.
"""
function assess_identifiability(
    ode::ModelingToolkit.ODESystem;
    measured_quantities = Array{ModelingToolkit.Equation}[],
    funcs_to_check = [],
    prob_threshold = 0.99,
    loglevel = Logging.Info,
)
    restart_logging(loglevel = loglevel)
    with_logger(_si_logger[]) do
        return _assess_identifiability(
            ode,
            measured_quantities = measured_quantities,
            funcs_to_check = funcs_to_check,
            prob_threshold = prob_threshold,
        )
    end
end

function _assess_identifiability(
    ode::ModelingToolkit.ODESystem;
    measured_quantities = Array{ModelingToolkit.Equation}[],
    funcs_to_check = [],
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

    result = _assess_identifiability(
        ode,
        funcs_to_check = funcs_to_check_,
        prob_threshold = prob_threshold,
    )
    nemo2mtk = Dict(funcs_to_check_ .=> funcs_to_check)
    out_dict = OrderedDict(nemo2mtk[param] => result[param] for param in funcs_to_check_)
    return out_dict
end

using PrecompileTools
include("precompile.jl")

end
