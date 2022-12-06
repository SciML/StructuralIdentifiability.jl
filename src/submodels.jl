


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





function find_submodels(ode)
    graph = construct_graph(ode)
    y = collect(keys(ode.y_equations))
    raw_models = traverse_outputs(graph, y)
    submodels = find_raw_submodels(raw_models, y, graph)
    unions = (search_add_unions(submodels))
    result = remove_empty(union(sort_all(unions)))
    return result
end