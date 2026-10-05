using DataStructures
using StructuralIdentifiability
using StructuralIdentifiability: cmp_states, cmp_lie_states, cmp_prefer_params, cmp_lie

comparators = OrderedDict(
    :default => cmp_prefer_params,
    :lie => cmp_lie,
    :states => cmp_states,
    :lie_states => cmp_lie_states,
)

include("../benchmarks.jl")

exclude_models = [:MAPK_5out, :MAPK_5out_bis, :MAPK_6out, :QWWC, :QY, :TumorHu, :TumorPillis, :LeukaemiaLeon2021, :cLV2, :Covid2, :Ovarian_follicle, :TranAlRadhawi, :NFkB]

for (label, model) in benchmarks
    if label in exclude_models
        continue
    end
    println("Processing $(model[:name])")

    for (name, comparator) in comparators
        println("With $name")
        res = reparametrize_global(model[:ode], cmp = comparator(model[:ode]))
        println("New ode: \n $(res[1]) \n")
        println("New vars: \n $(res[2]) \n")
        println("Relations: $(res[3])\n")
        println("=========")
    end
    println("\n===============================")
end
