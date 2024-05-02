using StructuralIdentifiability
using StructuralIdentifiability: initial_identifiable_functions, @ODEmodel, dennums_to_fractions

include("../benchmarks.jl")

for model in benchmarks
    if !get(model, :skip, false)
        (idfuncs_nostates, _) = initial_identifiable_functions(model[:ode], prob_threshold=0.9)
        push!(
             cases_simplification, 
             Dict(
                 :description => model[:name] * " io-coefficients",
                 :gens => idfuncs_nostates,
             )
        )
        (idfuncs_states, _) = initial_identifiable_functions(model[:ode], with_states = true, prob_threshold=0.9)
        push!(
            cases_simplification,
            Dict(
                :description => model[:name] * " simplified io-coeffs and Lie ders",
                :gens => idfuncs_states,
            )
        )
    end
end
