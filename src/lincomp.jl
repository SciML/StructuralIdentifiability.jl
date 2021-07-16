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
    linear_compartment_model(graph, inputs, outputs, leaks)

Input: defines a linear compartment model with nodes numbered from 1 to `n` by
- `graph` - and array of integer arrays representing the adjacency lists of the graph
- `inputs` - array of input nodes
- `outputs` - array of output nodes
- `leaks` - array of sink nodes

Output:
- the corresponding ODE system in the notation of https://doi.org/10.1007/s11538-015-0098-0
"""
function linear_compartment_model(
    graph::Vector{Vector{Int}},
    inputs::Vector{Int},
    outputs::Vector{Int},
    leaks::Vector{Int}
)
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

    R, vars = StructuralIdentifiability.Nemo.PolynomialRing(
        StructuralIdentifiability.Nemo.QQ, 
        vcat(x_vars_names, y_vars_names, u_vars_names, edges_vars_names)
    )
    x_vars = @view vars[1:n]
    x_equations = Dict{fmpq_mpoly, Union{fmpq_mpoly, Generic.Frac{fmpq_mpoly}}}(x => R(0) for x in x_vars)
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

    y_equations = Dict{fmpq_mpoly, Union{fmpq_mpoly, Generic.Frac{fmpq_mpoly}}}(str_to_var("y$i", R) => str_to_var("x$i", R) for i in outputs)

    return ODE{fmpq_mpoly}(x_equations, y_equations, Array{fmpq_mpoly}([str_to_var("u$i", R) for i in inputs]))
end

#------------------------------------------------------------------------------
