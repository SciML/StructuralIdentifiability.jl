"""
    saturate_outputs(ode, orders, max_deg = 5)

Adds extra outputs extracted from the Lie derivatives of the existing outputs up to the prescribed order.
Adds only the outputs which are algebraically independent over the existing outputs, parameters, and inputs
(and, thus, may be useful in the differential elimination process).

Inputs:
- `ode` - ode model
- `orders` - dictionary from outputs to the order of the Lie derivative to consider
- `max_deg` - the maximal degree (as a rational function) of possible new outputs to consider

Output: a new ode model with inputs inserted
"""
function saturate_outputs(
        ode::ODE{P},
        orders::Dict{P, Int},
        max_deg = 5,
    ) where {P <: MPolyRingElem}
    lie_ders = lie_derivatives_up_to(ode, Dict(ode.y_equations[y] => ord for (y, ord) in orders))
    @info "Degrees: $([(total_degree(numerator(f)), total_degree(denominator(f))) for f in lie_ders])"

    current_y = RationalFunctionField(vcat(collect(values(ode.y_equations)), ode.u_vars, ode.parameters))
    lie_ders = filter(f -> total_degree_frac(f) <= max_deg, lie_ders)
    sort!(lie_ders, lt = RationalFunctionFields.rational_function_cmp)
    for f in lie_ders
        if !first(RationalFunctionFields.check_algebraicity_modp(current_y, [f]))
            current_y = RationalFunctionField([generators(current_y)..., f])
        end
    end
    new_outputs = generators(current_y)[(length(ode.y_vars) + length(ode.u_vars) + length(ode.parameters) + 1):end]

    idx = 1
    old_y_names = map(var_to_str, ode.y_vars)
    new_ys = Dict{String, Generic.FracFieldElem{P}}()
    @info "New outputs $new_outputs"
    for new_y in new_outputs
        while "y_aux_$idx(t)" in old_y_names
            idx += 1
        end
        new_ys["y_aux_$idx(t)"] = new_y
        idx += 1
    end
    @info new_ys

    return add_outputs(ode, new_ys)
end

#------------------------------------------------------------------------------
"""
    propose_orders(ode)

For a given ode, proposes (using a heuristic) the orders for the outputs to consider in the saturation process.
The order is defined a the highest order of differentiation at which a new output may occur plus two.

Inputs:
- `ode` - ode model

Output: dictionary from outputs to orders
"""
function propose_orders(ode)
    x_with_u = filter(x -> length(intersect(vars(ode.x_equations[x]), ode.u_vars)) > 0, ode.x_vars)

    graph = Dict(x => Set(intersect(vars(ode.x_equations[x]), ode.x_vars)) for x in ode.x_vars)

    result = Dict(y => 0 for y in ode.y_vars)
    for y in ode.y_vars
        current_x = Set(intersect(vars(ode.y_equations[y]), ode.x_vars))
        reached_u = length(intersect(current_x, x_with_u))
        if reached_u > 0
            result[y] = 2
        end
        step = 1
        while true
            new_x = copy(current_x)
            for x in current_x
                new_x = union(new_x, graph[x])
            end
            if current_x == new_x
                break
            end
            current_x = new_x
            new_reached = length(intersect(current_x, x_with_u))
            if new_reached > reached_u
                result[y] = step + 2
            end
            reached_u = new_reached
            step += 1
        end
    end

    return result
end
