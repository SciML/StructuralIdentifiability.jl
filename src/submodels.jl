

"""
power_series_solution(ode, param_values, initial_conditions, input_values, prec)

Input:
- `ode` - an ode to solve
- `param_values` - parameter values, must be a dictionary mapping parameter to a value
- `initial_conditions` - initial conditions of `ode`, must be a dictionary mapping state variable to a value
- `input_values` - power series for the inpiuts presented as a dictionary variable => list of coefficients
- `prec` - the precision of solutions

Output: 
- computes a power series solution with precision prec presented as a dictionary variable => corresponding coordiante of the solution
"""

function construct_graph(ode) # construct a graph from given system of equations
    graph = Dict()
    for (x,f) in ode.x_equations
        graph[x] = vars(f)
    end
    for (y,f) in ode.y_equations
        graph[y] = vars(f)
    end
    return graph
end


function dfs(graph, start, visited) # dfs algo for traversing obtained graph
    push!(visited, start)
    if start in keys(graph)
        for node in graph[start]
            if !(node in visited)
                dfs(graph, node, visited)
            end 
        end
    end
    return visited
end



function traverse_outputs(graph, y)
    raw_models = Dict()
    for state in y
        model = dfs(graph, state, [])
        raw_models[state] = model
    end
    return raw_models
end


function find_raw_submodels(raw_models,Y, graph)
    result = []
    for y1 in Y
        for y2 in Y
            if y1 != y2
                if issubset(graph[y2], raw_models[y1]) && !(y2 in raw_models[y1])
                   push!(raw_models[y1], y2)
                end
            end
        end
    end
    for y in Y
        push!(result, raw_models[y])
    end
    return result
end


function sort_all(submodels)
    sorted = []
    for submodel in submodels
        sub = sort(submodel, by = string)
        push!(sorted, sub)
    end
    return sorted
end

function search_add_unions(submodels)
    result = [[]]
    for model in submodels
        for index in range(1,length(result))
            candidate_model = union(result[index], model)
            push!(result, candidate_model)
        end
    end
    return result
end
    

function remove_empty(submodels)
    if [] in submodels
        deleteat!(submodels, findfirst(x -> x == [], submodels))
    end
    return submodels
end





"""
    find_submodels(ode)

Input: 
- 'ode'

Output: 
- a tuple consisting of the power series solution and a dictionary of the form `(u, v) => power series`, where `u` is a state variable 
  `v` is a state or parameter, and the power series is the partial derivative of
  the function `u` w.r.t. `v` evaluated at the solution
"""

function submodel2ode_aux(ode, submodel)
    #print(ode)
    # Gleb: a more natural code would be to start with empty dicts and add the appropriate equations
    new_y = copy(ode.y_equations)
    new_x = copy(ode.x_equations)
    new_u = copy(ode.u_vars)
    for (x,f) in ode.x_equations
        if !(issubset(vars(x),submodel) && issubset(vars(f),submodel))
            delete!(new_x, x)
        end

    end
    
    for (y,f) in ode.y_equations
        if !(issubset(vars(y),submodel) && issubset(vars(f),submodel))
            delete!(new_y, y)
        end
    end
    # Gleb: there is a subtlety, so it is better to import `var_to_str` from SI.jl
    # and do `map(var_to_str, submodel)`
    sub = map(var_to_str, submodel)
    # Gleb: when incorporating into SI.jl, they should be referred as Nemo.PolynomialRing and Nemo.QQ
    S, _ = PolynomialRing(QQ, sub)
    fin_x = Dict(parent_ring_change(x, S) => parent_ring_change(f, S) for (x,f) in new_x)
    fin_y = Dict(parent_ring_change(y, S) => parent_ring_change(f, S) for (y,f) in new_y)

    # check which elem of submodel are in u_vars, change ring for them and add to fin_u
    # Gleb: it now remains to figure out the inputs
    # the easiest way would be to just go over the constructed equations and check which of the 
    # inputs are there
    return ODE{fmpq_mpoly}(fin_x, fin_y, Array{fmpq_mpoly, 1}())
end



"""
    power_series_solution(ode, param_values, initial_conditions, input_values, prec)

Input:
- `ode` - an ode to solve
- `param_values` - parameter values, must be a dictionary mapping parameter to a value
- `initial_conditions` - initial conditions of `ode`, must be a dictionary mapping state variable to a value
- `input_values` - power series for the inpiuts presented as a dictionary variable => list of coefficients
- `prec` - the precision of solutions

Output: 
- computes a power series solution with precision prec presented as a dictionary variable => corresponding coordiante of the solution
"""
# Gleb: this can be done in a single comprehension
# [ode_aux(ode2, submodel) for submodel in submodels]
function submodel2ode(ode2, submodels)
    return [submodel2ode_aux(ode2, submodel) for submodel in submodels]
end



"""
    power_series_solution(ode, param_values, initial_conditions, input_values, prec)

Input:
- `ode` - an ode to solve
- `param_values` - parameter values, must be a dictionary mapping parameter to a value
- `initial_conditions` - initial conditions of `ode`, must be a dictionary mapping state variable to a value
- `input_values` - power series for the inpiuts presented as a dictionary variable => list of coefficients
- `prec` - the precision of solutions

Output: 
- computes a power series solution with precision prec presented as a dictionary variable => corresponding coordiante of the solution
"""


function visualize_graph(ode)
    graph = construct_graph(ode)
    y = ode.y_vars
    x = ode.x_vars
    u = ode.u_vars
    states = vcat(y,x,u)
    member_y = [1 for i in y]
    member_x = [2 for i in x]
    member_u = [3 for i in u]
    member = vcat(member_y, member_x, member_u)
    nodecolor = [colorant"lightblue", colorant"lightgreen", colorant"lightpink"]
    nodefillc = nodecolor[member]
    dict = Dict()
    for (i,el) in enumerate(states)
        dict[el] = i
    end

    dict

    g = SimpleDiGraph(length(states))

    for (x,f) in graph
        for node in f
            if (node != x && (node in states))
                add_edge!(g, dict[x], dict[node])
            end
        end
    end

    gplot(g, layout=circular_layout,nodefillc=nodefillc, nodelabel=states, arrowlengthfrac =0.05, EDGELINEWIDTH = 0.2, edgestrokec = colorant"grey", NODELABELSIZE =3)
end

"""
    power_series_solution(ode, param_values, initial_conditions, input_values, prec)

Input:
- `ode` - an ode to solve
- `param_values` - parameter values, must be a dictionary mapping parameter to a value
- `initial_conditions` - initial conditions of `ode`, must be a dictionary mapping state variable to a value
- `input_values` - power series for the inpiuts presented as a dictionary variable => list of coefficients
- `prec` - the precision of solutions

Output: 
- computes a power series solution with precision prec presented as a dictionary variable => corresponding coordiante of the solution
"""

function find_submodels(ode)
    graph = construct_graph(ode)
    y = collect(keys(ode.y_equations))
    raw_models = traverse_outputs(graph, y)
    submodels = find_raw_submodels(raw_models, y, graph)
    unions = (search_add_unions(submodels))
    result = remove_empty(union(sort_all(unions)))
    back2ode = submodel2ode(ode, result)
    return back2ode
end