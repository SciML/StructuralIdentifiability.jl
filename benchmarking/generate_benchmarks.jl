using Printf

using StructuralIdentifiability
using StructuralIdentifiability: ODE, print_for_DAISY, print_for_COMBOS, print_for_maple

include("benchmarks.jl")

formats = [
    Dict(
        :name => "SIAN",
        :function => ode -> print_for_maple(ode, :SIAN),
        :extension => ".mpl"
    ),
    Dict(
        :name => "RosenfeldGroebner",
        :function => ode -> print_for_maple(ode, :DifferentialAlgebra),
        :extension => ".mpl"
    ),
    Dict(
        :name => "DifferentialThomas",
        :function => ode -> print_for_maple(ode, :DifferentialThomas),
        :extension => ".mpl"
    ),
    Dict(
        :name => "DAISY",
        :function => print_for_DAISY,
        :extension => ".txt"
    ),
    Dict(
        :name => "COMBOS",
        :function => print_for_COMBOS,
        :extension => ".txt"
    )
]

for frmt in formats
    mkpath(frmt[:name])
    cd(frmt[:name])
    for bnchmrk in benchmarks
        as_text = frmt[:function](bnchmrk[:ode])
        fname = replace(bnchmrk[:name], " " => "-") * frmt[:extension]
        open(fname, "w") do io
            write(io, as_text)
        end
    end
    cd("..")
end
