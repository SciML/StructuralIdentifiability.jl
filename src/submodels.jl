

"""
    construct_graph(ode)

This function creates a dictionary with all the dependencies in terms of variables expressed in the system of ODEs given.

Input:
- `ode` - an ODEs system 

Output: 
- Dictionary where each key is a variable and each value is a list of variables on which the key depends
"""

function construct_graph(ode) 
    graph = Dict()
    for (x,f) in ode.x_equations
        graph[x] = vars(f)
    end
    for (y,f) in ode.y_equations
        graph[y] = vars(f)
    end
    return graph
end


function dfs(graph, start, visited) 
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


function traverse_outputs(graph, ys)
    raw_models = Dict()
    for y in ys
        model = dfs(graph, y, [])
        raw_models[y] = model
    end
    return raw_models
end


function find_raw_submodels(raw_models, ys, graph)
    result = []
    for y1 in ys
        for y2 in ys
            if y1 != y2
                if issubset(graph[y2], raw_models[y1]) && !(y2 in raw_models[y1])
                   push!(raw_models[y1], y2)
                end
            end
        end
    end
    for y in ys
        push!(result, raw_models[y])
    end
    return result
end

# Gleb: why do we need this ?
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
        for index in 1:length(result)
            push!(result, union(result[index], model))
        end
    end
    return result
end
    

function ode_aux(ode, submodel)
    new_y = copy(ode.y_equations)
    new_x = copy(ode.x_equations)
    new_u = Array{fmpq_mpoly, 1}([u for u in ode.u_vars if u in submodel])
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

    sub_str = map(var_to_str, submodel)
    S, _ = Nemo.PolynomialRing(Nemo.QQ, sub_str)
    fin_x = Dict(parent_ring_change(x, S) => parent_ring_change(f, S) for (x,f) in new_x)
    fin_y = Dict(parent_ring_change(y, S) => parent_ring_change(f, S) for (y,f) in new_y)
    fin_u = [parent_ring_change(u, S) for u in new_u]

    return ODE{fmpq_mpoly}(fin_x, fin_y, fin_u)
end


"""
    submodel2ode(ode, submodels)
This function converts valid submodels which are initially represented as lists of involved variables into ODE objects 
using the function `ode_aux` where one constructs a new system of ODEs for each valid submodel separately. 

Input:
- `ode` - initial ODE system
- `submodels` - list of valid submodels


Output: 
- A list of ODE objects, each corresponding to a certain valid submodel
"""
# Gleb: not sure if this list comprehension has to be a separate function
function submodel2ode(ode, submodels)
    return [ode_aux(ode, submodel) for submodel in submodels]
end



#"""
#    visualize_ode(ode)
#        
#This function produces a plot of the system of ODEs given in terms of its graph representatinon. 
#Each node is either a state, an observation or a constant and each directed edge e=(i,j) indicates 
#that the law of the variable i depends on the law of the variable j. 
#
#Input:
#- `ode` - an ODEs system
#
#Output: 
#- A temporary .html file (that is opened automatically) where the graph is displayed 
#"""
#function visualize_ode(ode)
#    graph = construct_graph(ode)
#    y = ode.y_vars
#    x = ode.x_vars
#    u = ode.u_vars
#    states = vcat(y,x,u)
#    member_y = [1 for i in y]
#    member_x = [2 for i in x]
#    member_u = [3 for i in u]
#    member = vcat(member_y, member_x, member_u)
#    nodecolor = [colorant"lightblue", colorant"lightgreen", colorant"lightpink"]
#    nodefillc = nodecolor[member]
#    ind = Dict()
#    for (i,el) in enumerate(states)
#        ind[el] = i
#    end
#
#    g = SimpleDiGraph(length(states))
#
#    for (x,f) in graph
#        for node in f
#            if (node != x && (node in states))
#                add_edge!(g, ind[x], ind[node])
#            end
#        end
#    end
#    gplothtml(g, layout=circular_layout,nodefillc=nodefillc, nodelabel=states, arrowlengthfrac =0.05, EDGELINEWIDTH = 0.2, edgestrokec = colorant"grey", NODELABELSIZE =3)
#end


# TODO: add an example to this docstring
"""
    find_submodels(ode)
 
The function calculates and returns all valid submodels given a system of ODEs.

Input:
- `ode` - an ODEs system to be studied

Output: 
- A list of submodels represented as `ode` objects
"""
function find_submodels(ode)
    graph = construct_graph(ode)
    y = collect(keys(ode.y_equations))
    raw_models = traverse_outputs(graph, y)
    submodels = find_raw_submodels(raw_models, y, graph)
    unions = search_add_unions(submodels)
    result = [m for m in union(sort_all(unions)) if length(m) > 0]
    back2ode = submodel2ode(ode, result)
    return back2ode
end
