# Copyright (c) 2021, R. Dong, C. Goodbreak, H. Harrington, G. Pogudin
# Copyright (c) 2020, A. Ovchinnikov, A. Pillay, G. Pogudin, T. Scanlon

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#------------------------------------------------------------------------------

"""
    linear_compartment_model(graph; inputs = [], outputs = [], leaks = [])

Input: defines a linear compartment model with nodes numbered from 1 to `n` by
- `graph` - and array of integer arrays representing the adjacency lists of the graph
- `inputs` - array of input nodes
- `outputs` - array of output nodes
- `leaks` - array of sink nodes

Output:
- the corresponding ODE system in a standard notation  (as, e.g., in [this paper](https://doi.org/10.1007/s11538-015-0098-0))

Example: Consider a bidirected cycle with four nodes. Its adjacency list can be written as follows:
```
[ [2, 4], [1, 3], [2, 4], [1, 3] ]
```
In the list above, the `i`-th element is a list of vertices to which there exists an edge
from the vertex `i`. Now we can create a linear compartment model over this graph with
the output at vertex 1, input at vertex 2, and leaks at vertices 3 and 4 as follows:
```jldoctest
julia> ode = linear_compartment_model([[2, 4], [1, 3], [2, 4], [1, 3]], outputs = [1], inputs = [2], leaks = [2, 3])
x1' = -x1*a_2_1 - x1*a_4_1 + x2*a_1_2 + x4*a_1_4
x3' = x2*a_3_2 - x3*a_2_3 - x3*a_4_3 - x3*a_0_3 + x4*a_3_4
x2' = x1*a_2_1 - x2*a_1_2 - x2*a_3_2 - x2*a_0_2 + x3*a_2_3 + u2
x4' = x1*a_4_1 + x3*a_4_3 - x4*a_1_4 - x4*a_3_4
y1 = x1
```
"""
function linear_compartment_model(graph; inputs = [], outputs = [], leaks = [])
    n = length(graph)
    x_vars_names = ["x$i" for i in 1:n]
    y_vars_names = ["y$i" for i in outputs]
    u_vars_names = ["u$i" for i in inputs]
    edges_vars_names = Array{String, 1}()
    for i in 1:n
        for j in graph[i]
            push!(edges_vars_names, "a_$(j)_$(i)")
        end
    end
    for s in leaks
        push!(edges_vars_names, "a_0_$(s)")
    end

    R, vars = StructuralIdentifiability.Nemo.polynomial_ring(
        StructuralIdentifiability.Nemo.QQ,
        vcat(x_vars_names, y_vars_names, u_vars_names, edges_vars_names),
    )
    x_vars = @view vars[1:n]
    x_equations =
        Dict{QQMPolyRingElem, ExtendedFraction{QQMPolyRingElem}}(x => R(0) for x in x_vars)
    for i in 1:n
        for j in graph[i]
            rate = str_to_var("a_$(j)_$(i)", R)
            x_equations[x_vars[j]] += x_vars[i] * rate
            x_equations[x_vars[i]] -= x_vars[i] * rate
        end
        if i in leaks
            rate = str_to_var("a_0_$(i)", R)
            x_equations[x_vars[i]] += -x_vars[i] * rate
        end
        if i in inputs
            x_equations[x_vars[i]] += str_to_var("u$i", R)
        end
    end

    y_equations = Dict{QQMPolyRingElem, ExtendedFraction{QQMPolyRingElem}}(
        str_to_var("y$i", R) => str_to_var("x$i", R) for i in outputs
    )

    return ODE{QQMPolyRingElem}(
        x_equations,
        y_equations,
        Array{QQMPolyRingElem}([str_to_var("u$i", R) for i in inputs]),
    )
end

#------------------------------------------------------------------------------
