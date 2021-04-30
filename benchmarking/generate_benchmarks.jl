using Printf

using StructuralIdentifiability
using StructuralIdentifiability: ODE, print_for_DAISY, print_for_SIAN

include("benchmarks.jl")

formats = [
    Dict(
        :name => "SIAN",
        :function => print_for_SIAN,
        :extension => ".mpl"
    ),
    Dict(
        :name => "DAISY",
        :function => print_for_DAISY,
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
