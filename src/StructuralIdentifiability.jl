module StructuralIdentifiability

# General purpose packages
using Base.Iterators
using Combinatorics
using DataStructures
using Dates
using IterTools
using LinearAlgebra
using Logging
using MacroTools
using Primes
using Random

# Algebra packages
using AbstractAlgebra
using Nemo
using Groebner
using ParamPunPam
using ParamPunPam: reduce_mod_p!, specialize_mod_p

using ModelingToolkit

# defining a model
export ODE, @ODEmodel, preprocess_ode

# assessing identifiability
export assess_local_identifiability, assess_global_identifiability, assess_identifiability

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

# would be great to merge with the Julia logger
const _runtime_logger = Dict(:id_calls_to_gb => 0, :id_groebner_time => 0.0)

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
include("RationalFunctionFields/IdealMQS.jl")
include("RationalFunctionFields/RationalFunctionField.jl")
include("global_identifiability.jl")
include("identifiable_functions.jl")
include("lincomp.jl")
include("pb_representation.jl")
include("submodels.jl")
include("discrete.jl")

# TODO: print equations in ODEmodel in the order they were given in the macro
# TODO: handle finding identifiabile functions in case there are no parameters

"""
    assess_identifiability(ode::ODE{P}, p::Float64=0.99) where P <: MPolyElem{fmpq}

Input:
- `ode` - the ODE model
- `p` - probability of correctness.

Assesses identifiability (both local and global) of a given ODE model (parameters detected automatically). The result is guaranteed to be correct with the probability
at least `p`.
"""
function assess_identifiability(ode::ODE{P}, p::Float64 = 0.99) where {P <: MPolyElem{fmpq}}
    result = assess_identifiability(ode, vcat(ode.parameters, ode.x_vars), p)
    return result #Dict(param => res for (param, res) in zip(ode.parameters, result_list))
end

"""
    assess_identifiability(ode, [funcs_to_check, p=0.99])

Input:
- `ode` - the ODE model
- `p` - probability of correctness.
    
Assesses identifiability of a given ODE model. The result is guaranteed to be correct with the probability
at least `p`.

If `funcs_to_check` are given, then the function will assess the identifiability of the provided functions
and return a list of the same length with each element being one of `:nonidentifiable`, `:locally`, `:globally`.

If `funcs_to_check` are not given, the function will assess identifiability of the parameters, and the result will
be a dictionary from the parameters to their identifiability properties (again, one of `:nonidentifiable`, `:locally`, `:globally`).
"""
function assess_identifiability(
    ode::ODE{P},
    funcs_to_check::Array{<:RingElem, 1},
    p::Float64 = 0.99,
) where {P <: MPolyElem{fmpq}}
    p_glob = 1 - (1 - p) * 0.9
    p_loc = 1 - (1 - p) * 0.1

    @info "Assessing local identifiability"
    trbasis = Array{fmpq_mpoly, 1}()
    runtime = @elapsed local_result =
        assess_local_identifiability(ode, funcs_to_check, p_loc, :SE, trbasis)
    @info "Local identifiability assessed in $runtime seconds"
    @debug "Trasncendence basis to be specialized is $trbasis"
    _runtime_logger[:loc_time] = runtime

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

    result = Dict{Any, Symbol}()
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
    assess_identifiability(ode::ModelingToolkit.ODESystem; measured_quantities=Array{ModelingToolkit.Equation}[], funcs_to_check=[], p = 0.99)

Input:
- `ode` - the ModelingToolkit.ODESystem object that defines the model
- `measured_quantities` - the output functions of the model
- `funcs_to_check` - functions of parameters for which to check the identifiability
- `p` - probability of correctness.

Assesses identifiability (both local and global) of a given ODE model (parameters detected automatically). The result is guaranteed to be correct with the probability
at least `p`.
"""
function assess_identifiability(
    ode::ModelingToolkit.ODESystem;
    measured_quantities = Array{ModelingToolkit.Equation}[],
    funcs_to_check = [],
    p = 0.99,
)
    if isempty(measured_quantities)
        measured_quantities = get_measured_quantities(ode)
    end
    if length(funcs_to_check) == 0
        funcs_to_check = ModelingToolkit.parameters(ode)
    end
    ode, conversion = preprocess_ode(ode, measured_quantities)
    out_dict = Dict{Num, Symbol}()
    funcs_to_check_ = [eval_at_nemo(each, conversion) for each in funcs_to_check]
    result = assess_identifiability(ode, funcs_to_check_, p)
    nemo2mtk = Dict(funcs_to_check_ .=> funcs_to_check)
    out_dict = Dict(nemo2mtk[param] => result[param] for param in funcs_to_check_)
    return out_dict
end

using PrecompileTools
include("precompile.jl")

end
