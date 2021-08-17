module StructuralIdentifiability

# General purpose packages
using Base.Iterators
using Dates
using IterTools
using LinearAlgebra
using Logging
using MacroTools
using Primes
using DataStructures

# Algebra packages
using AbstractAlgebra
using Nemo
using GroebnerBasis
using Singular

# For testing (TODO: move to the test-specific dependencies)
using Test
using TestSetExtensions
using IterTools
using MacroTools
using ModelingToolkit

# defining a model
export ODE, @ODEmodel, PreprocessODE

# assessing identifiability
export assess_local_identifiability, assess_global_identifiability, assess_identifiability

# auxuliary function
export set_parameter_values

# extra functionality
export find_ioequations, find_identifiable_functions

# exporting to other formats
export print_for_maple, print_for_DAISY

# function for creating linear compartment models
export linear_compartment_model

# would be great to merge with the Julia logger
_runtime_logger = Dict()

include("util.jl")
include("power_series_utils.jl")
include("ODE.jl")
include("ODEexport.jl")
include("local_identifiability.jl")
include("wronskian.jl")
include("elimination.jl")
include("primality_check.jl")
include("io_equation.jl")
include("global_identifiability.jl")
include("lincomp.jl")

"""
    assess_identifiability(ode::ODE{P}, p::Float64=0.99) where P <: MPolyElem{fmpq}
Input:
- `ode` - the ODE model
- `p` - probability of correctness.

Assesses identifiability (both local and global) of a given ODE model (parameters detected automatically). The result is guaranteed to be correct with the probability
at least `p`.

"""
function assess_identifiability(ode::ODE{P}, p::Float64=0.99) where P <: MPolyElem{fmpq}
    result_list = assess_identifiability(ode, ode.parameters, p)
    return Dict(param => res for (param, res) in zip(ode.parameters, result_list))
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
function assess_identifiability(ode::ODE{P}, funcs_to_check::Array{<: RingElem,1}, p::Float64=0.99) where P <: MPolyElem{fmpq} 
    p_glob = 1 - (1 - p) * 0.9
    p_loc = 1 - (1 - p) * 0.1

    @info "Assessing local identifiability"
    runtime = @elapsed local_result, bound = assess_local_identifiability(ode, funcs_to_check, p_loc, :ME)
    @info "Local identifiability assessed in $runtime seconds"
    _runtime_logger[:loc_time] = runtime

    if bound > 1
        @warn "For this model single-experiment identifiable functions are not the same as multi-experiment identifiable. The analysis performed by this program is right now for multi-experiment only. If you still  would like to assess single-experiment identifiability, we recommend using SIAN (https://github.com/alexeyovchinnikov/SIAN-Julia)"
        @debug "Bound: $bound"
    end

    locally_identifiable = Array{Any,1}()
    for (loc, f) in zip(local_result, funcs_to_check)
        if loc
            push!(locally_identifiable, f)
        end
    end

    @info "Assessing global identifiability"
    runtime = @elapsed global_result = assess_global_identifiability(ode, locally_identifiable, p_glob)
    @info "Global identifiability assessed in $runtime seconds" 
    _runtime_logger[:glob_time] = runtime

    result = Array{Symbol,1}()
    glob_ind = 1
    for i in 1:length(funcs_to_check)
        if !local_result[i]
            push!(result, :nonidentifiable)
        else
            if global_result[glob_ind]
                push!(result, :globally)
            else
                push!(result, :locally)
            end
            glob_ind += 1
        end
    end

    return result
end

function assess_identifiability(ode::ModelingToolkit.ODESystem, inputs, funcs_to_check, p::Float64=0.99)
    diff_eqs = equations(ode)
    params = ModelingToolkit.parameters(ode)
    state_vars = ModelingToolkit.states(ode)
    y_functions = [each.lhs for each in ModelingToolkit.observed(ode)]
    output_eqs = ModelingToolkit.observed(ode) 
    ode, syms, gens_ = PreprocessODE(diff_eqs, output_eqs, state_vars, y_functions, inputs, params)
    if length(funcs_to_check) > 0
        funcs_to_check = [substitute(x, syms .=> gens_) for x in funcs_to_check]
        return assess_identifiability(ode, funcs_to_check, p)
    else
        return assess_identifiability(ode, p)
    end
end

end
