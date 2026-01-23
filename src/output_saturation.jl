function saturate_outputs(
    ode::ODE{P},
    orders::Dict{P, Int},
) where {P <: MPolyRingElem}
    lie_ders = lie_derivatives_up_to(ode, Dict(ode.y_equations[y] => ord for (y, ord) in orders))
    @info lie_ders
    @info "Degrees: $([(total_degree(numerator(f)), total_degree(denominator(f))) for f in lie_ders])"

    current_y = RationalFunctionField(collect(values(ode.y_equations)))
    lie_ders = filter(f -> total_degree_frac(f) <= 5, lie_ders)
    sort!(lie_ders, lt = RationalFunctionFields.rational_function_cmp)
    for f in lie_ders
        if !first(RationalFunctionFields.field_contains_mod_p(current_y, [f]))
            current_y = RationalFunctionField([generators(current_y)..., f])
        end
    end
    new_outputs = generators(current_y)[length(ode.y_vars) + 1 : end]

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

function propose_orders(ode)
    x_with_u = filter(x -> length(intersect(vars(ode.x_equations[x]), ode.u_vars)) > 0, ode.x_vars)

    graph = Dict(x => Set(intersect(vars(ode.x_equations[x]), ode.x_vars)) for x in ode.x_vars)

    result = Dict(y => 0 for y in ode.y_vars)
    for y in ode.y_vars
        current_x = Set(intersect(vars(ode.y_equations[y]), ode.x_vars))
        reached_u = length(intersect(current_x, x_with_u))
        if reached_u > 0
            result[y] = 1
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
                result[y] = step + 1
            end
            reached_u = new_reached
            step += 1
        end
    end

    return result
end
