using StructuralIdentifiability, RationalFunctionFields
using StructuralIdentifiability: lie_derivatives_up_to

function simplify_lie(ode, order; simplify = :weak)
    orders = Dict(
        ode.y_equations[y] => order
        for y in ode.y_vars
    )
    generators = lie_derivatives_up_to(ode, orders)
    simpler = simplified_generating_set(
        RationalFunctionField(generators);
        simplify = simplify,
        prob_threshold = 0.99,
    )
    generators, simpler
end
