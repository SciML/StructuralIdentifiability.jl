using HomotopyContinuation

"""
    Everything for locally identifiable, polynomial, single-output, no-input models for the moment
"""

"""
    Computes enough prolongations by Lie derivatives to use SIAN
"""
function get_poly_system_Lie(ode::ODE{P}) where P <: MPolyElem{<: FieldElem}
    prolong_count = length(ode.x_vars) + length(ode.parameters)

    # creating differential ring
    R = DiffPolyRing(
        Nemo.QQ, 
        map(var_to_str, vcat(ode.x_vars, ode.y_vars)), 
        map(var_to_str, ode.parameters),
        vcat([0 for x in ode.x_vars], [prolong_count for y in ode.y_vars])
    )
    set_custom_derivations!(R, ode.x_equations)

    # equations forming square subsystem
    result_core = Array{P, 1}()
    # extra equations
    result_extra = Array{P, 1}()

    for y in ode.y_vars
        eq = y - ode.y_equations[y]
        push!(result_core, to_diffpoly(eq, R))
        for _ in 1:(prolong_count - 1)
            push!(result_core, diff(result_core[end], R))
        end
        push!(result_extra, diff(result_core[end], R))
    end

    return Dict(
        :diff_ring => R,
        :square_system => result_core, 
        :extra_equations => result_extra
    )
end

#------------------------------------------------------------------------------

"""
    Converts an ODE into a digraph sending each variable in lhs to the variables on the rhs
"""
function ode_to_graph(ode::ODE{P}) where P <: MPolyElem{<: FieldElem}
    graph = Dict()

    for (z, eq) in merge(ode.y_equations, ode.x_equations)
        graph[var_to_str(z)] = [var_to_str(v) for v in vars(eq) if !(v in ode.parameters)]
    end
    for u in ode.u_vars
        graph[var_to_str(u)] = Array{String, 1}()
    end

    return graph
end

#------------------------------------------------------------------------------

"""
    Finds the prolongation order for all the variables based on the desired orders of the
    outputs via DFS of the corresponding graph
"""
function get_prolongation_orders(graph, out_orders)
    result = copy(out_orders)
    stack = collect(keys(result))
    for v in keys(graph)
        if !(v in keys(result))
            result[v] = -1
        end
    end

    while length(stack) > 0
        v = pop!(stack)
        offset = (v in keys(out_orders)) ? 0 : -1
        for u in graph[v]
            if result[u] < result[v] + offset
                result[u] = result[v] + offset
                push!(stack, u)
            end
        end
    end

    return result
end

#------------------------------------------------------------------------------

"""
    Computes enough prolongations by the lazy approach to use SIAN
"""
function get_poly_system_lazy(ode::ODE{P}) where P <: MPolyElem{<: FieldElem}
    prolong_count = length(ode.x_vars) + length(ode.parameters)

    graph = ode_to_graph(ode)
    prolongation_total = get_prolongation_orders(
        graph, 
        Dict(var_to_str(y) => prolong_count for y in ode.y_vars)
    )
    prolongation_square = get_prolongation_orders(
        graph,
        Dict(var_to_str(y) => prolong_count - 1 for y in ode.y_vars)
    )

    # creating differential ring
    diff_vars = map(var_to_str, vcat(ode.x_vars, ode.y_vars, ode.u_vars))
    orders = [prolongation_total[v] for v in diff_vars]
    R = DiffPolyRing(
        Nemo.QQ, 
        diff_vars,
        map(var_to_str, ode.parameters),
        orders
    )

    # equations forming square subsystem
    result_core = Array{P, 1}()
    # extra equations
    result_extra = Array{P, 1}()

    for y in ode.y_vars
        eq = y - ode.y_equations[y]
        push!(result_core, to_diffpoly(eq, R))
        for _ in 1:prolongation_square[var_to_str(y)]
            push!(result_core, diff(result_core[end], R))
        end
        push!(result_extra, diff(result_core[end], R))
    end

    for x in ode.x_vars
        if prolongation_total[var_to_str(x)] <= 0
            break
        end
        eq = str_to_var(diffvar(var_to_str(x), 1), R.ring) - to_diffpoly(ode.x_equations[x], R)
        push!(result_core, eq)
        for _ in 2:prolongation_total[var_to_str(x)]
            push!(result_core, diff(result_core[end], R))
        end
        #append!(result_extra, prolongations[(prolongation_square[var_to_str(x)] + 1):end])
    end

    return Dict(
        :diff_ring => R,
        :square_system => result_core,
        :extra_equations => result_extra
    )
end

#------------------------------------------------------------------------------

function oscar_to_hc_poly(poly::MPolyElem, binding::Dict{<: MPolyElem, HomotopyContinuation.Variable})
    result = HomotopyContinuation.Expression(0)
    for (m, coef) in zip(exponent_vectors(poly), coeffs(poly))
        m_new = prod([binding[gens(parent(poly))[i]]^e for (i, e) in enumerate(m)])
        result += Float64(coef) * m_new
    end
    return result
end

#------------------------------------------------------------------------------

function numerical_identifiability(ode::ODE{P}, method=:lazy, solver=:monodromy) where P <: MPolyElem{<: FieldElem}
    # prolonging the system
    equations = undef
    if method == :lazy
        equations = get_poly_system_lazy(ode)
    elseif method == :Lie
        equations = get_poly_system_Lie(ode)
    else
        throw(Base.ArgumentError("Unknown method $method"))
    end
    R = equations[:diff_ring]

    # generating a solution
    param_vals = Dict(p => Nemo.QQ(rand(1:10) // 10) for p in ode.parameters)
    ic_vals = Dict(x => Nemo.QQ(rand(1:10) // 10) for x in ode.x_vars)
    ps_sol = power_series_solution(
        ode, 
        param_vals, 
        ic_vals, 
        Dict{P, Array{fmpq, 1}}(), length(ode.x_vars) + length(ode.parameters) + 2
    )
    ps_sol = Dict(var_to_str(v) => sol for (v, sol) in ps_sol)
    point = produce_point(ps_sol, R)

    y_varnames = map(var_to_str, ode.y_vars)
    params = Array{HomotopyContinuation.Variable, 1}()
    nonparams = collect(map(p -> Variable(var_to_str(p)), R.parameters))
    binding = Dict(p => Variable(var_to_str(p)) for p in R.parameters)

    for (i, v) in enumerate(R.diff_var_names)
        for ord in 0:R.max_orders[i]
            binding[str_to_var(diffvar(v, ord), R.ring)] = Variable(v, ord)
            if v in y_varnames 
                push!(params, Variable(v, ord))
            else
                push!(nonparams, Variable(v, ord))
            end
        end
    end
    sort!(params)
    sort!(nonparams)

    hc_polys = map(eq -> oscar_to_hc_poly(eq, binding), equations[:square_system])
    hc_polys_ext = map(eq -> oscar_to_hc_poly(eq, binding), vcat(equations[:square_system], equations[:extra_equations]))
    inv_binding = Dict(b => a for (a, b) in binding)
    point_param = [Float64(point[inv_binding[p]]) for p in params]
    point_nonparam = [Float64(point[inv_binding[np]]) for np in nonparams]

    sol_raw = undef
    if solver == :monodromy
        system = HomotopyContinuation.System(hc_polys; parameters=params)
        sol_raw = HomotopyContinuation.monodromy_solve(system, point_nonparam, point_param)
    elseif solver == :homotopy
        system = HomotopyContinuation.System([HomotopyContinuation.evaluate(eq, params => point_param) for eq in hc_polys])
        sol_raw = HomotopyContinuation.solve(System)
    else
        @info "No solver $solver"
    end

    @info sol_raw

    final_result = []
    for solution in solutions(sol_raw)
        to_add = true
        for poly in hc_polys_ext
            if abs(HomotopyContinuation.evaluate(poly, vcat(nonparams, params) => vcat(solution, point_param))) > 10^(-6)
                to_add = false
                println(HomotopyContinuation.evaluate(poly, vcat(nonparams, params) => vcat(solution, point_param)))
            end
        end
        if to_add
            push!(final_result, solution)
        end
    end

    @info "Number of solutions: $(length(final_result))"
    @info final_result

    return length(final_result) == 1
end
