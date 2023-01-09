# ------------------------------------------------------------------------------

"""
    construct_graph(ode)

This function creates a dictionary with all the dependencies in terms of variables expressed in the system of ODEs given.

Input:
- `ode` - an ODEs system 

Output: 
- Dictionary where each key is a variable and each value is a list of variables on which the key depends
"""

function construct_graph(ode::ODE{P}) where P <: MPolyElem
    graph = Dict{fmpq_mpoly, Set{fmpq_mpoly}}()
    for (x,f) in ode.x_equations
        temp = unpack_fraction(f)
        graph[x] = Set{fmpq_mpoly}(union(vars(temp[1]), vars(temp[2])))
    end
    for (y,f) in ode.y_equations
        temp = unpack_fraction(f)
        graph[y] = Set{fmpq_mpoly}(union(vars(temp[1]), vars(temp[2])))
    end
    
    return graph
end

# ------------------------------------------------------------------------------

function dfs(graph::Dict{fmpq_mpoly, Set{fmpq_mpoly}}, start::fmpq_mpoly, visited::Set{fmpq_mpoly})
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

# ------------------------------------------------------------------------------

function traverse_outputs(graph::Dict{fmpq_mpoly, Set{fmpq_mpoly}}, ys::Array{fmpq_mpoly, 1})
    raw_models = Dict{fmpq_mpoly, Set{fmpq_mpoly}}()
    for y in ys
        model = dfs(graph, y, Set{fmpq_mpoly}())
        raw_models[y] = model
    end
    return raw_models
end

# ------------------------------------------------------------------------------

function saturate_ys(unions::Array{Set{fmpq_mpoly}, 1},
                     Y::Array{fmpq_mpoly,1}, 
                     graph::Dict{fmpq_mpoly, Set{fmpq_mpoly}}, 
                     X::Array{fmpq_mpoly,1})
    for element in unions
        for y in Y
            states = [x for x in graph[y] if x in X]
            if issubset(states, element) && !(y in element)
                push!(element, y)
            end
        end
    end
end

# ------------------------------------------------------------------------------

function search_add_unions(submodels::Array{Set{fmpq_mpoly}, 1})
    result = Array{Set{fmpq_mpoly}, 1}([Set{fmpq_mpoly}()])
    for model in submodels
        for index in 1:length(result)
            push!(result, union(result[index], model))
        end
    end
    return result
end

# ------------------------------------------------------------------------------
   
function filter_max(ode::ODE{P}, submodels::Array{Set{fmpq_mpoly}, 1}) where P <: MPolyElem
    n = length(ode.x_vars)
    new_sub = Array{Set{fmpq_mpoly}, 1}([])
    for submodel in submodels
        list_x = [x for x in submodel if x in ode.x_vars]
        if !(length(list_x) == n)
            push!(new_sub, submodel)
        end
    end
    return new_sub
end
    
# ------------------------------------------------------------------------------

function ode_aux(ode::ODE{P}, submodel::Set{fmpq_mpoly}) where P <: MPolyElem
    new_y = copy(ode.y_equations)
    new_x = copy(ode.x_equations)
    new_u = Array{fmpq_mpoly, 1}([u for u in ode.u_vars if u in submodel])
    for (x,f) in ode.x_equations
        if !(issubset(vars(x), submodel) && issubset(vars(f), submodel))
            delete!(new_x, x)
        end
    end
    
    for (y,f) in ode.y_equations
        if !(issubset(vars(y), submodel) && issubset(vars(f), submodel))
            delete!(new_y, y)
        end
    end

    sub_str = map(var_to_str, collect(submodel))
    S, _ = Nemo.PolynomialRing(Nemo.QQ, sub_str)
    fin_x = Dict(parent_ring_change(x, S) => parent_ring_change(f, S) for (x,f) in new_x)
    fin_y = Dict(parent_ring_change(y, S) => parent_ring_change(f, S) for (y,f) in new_y)
    fin_u = [parent_ring_change(u, S) for u in new_u]

    return ODE{fmpq_mpoly}(fin_x, fin_y, fin_u)
end

# ------------------------------------------------------------------------------

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

function submodel2ode(ode::ODE{P}, submodels::Array{Set{fmpq_mpoly}, 1}) where P <: MPolyElem
    return [ode_aux(ode, submodel) for submodel in submodels]
end

# ------------------------------------------------------------------------------

"""
    find_submodels(ode)
 
The function calculates and returns all valid submodels given a system of ODEs.

Input:
- `ode` - an ODEs system to be studied

Output: 
- A list of submodels represented as `ode` objects
Example:
```
>ode = @ODEmodel(x1'(t) = x1(t)^2, 
                 x2'(t) = x1(t) * x2(t), 
                 y1(t) = x1(t), 
                 y2(t) = x2(t))
>find_submodels(ode)
    ODE{fmpq_mpoly}[
        
        x1'(t) = a(t)*x2(t)^2 + x1(t)
        y1(t) = x1(t)
    ]
```
"""
function find_submodels(ode::ODE{P}) where P <: MPolyElem
    
    graph = construct_graph(ode)
    ys = ode.y_vars
    xs = ode.x_vars
    raw_models = traverse_outputs(graph, ys)
    input_unions = [raw_models[y] for y in ys]
    unions = (search_add_unions(input_unions))
    saturate_ys(unions, ys, graph,xs)
    result = filter_max(ode, union(unions)[2:end])
    back2ode = submodel2ode(ode, result)
    return back2ode
end
