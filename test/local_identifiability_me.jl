# Copyright (c) 2021, R. Dong, C. Goodbarke, H. Harrington, G. Pogudin
# Copyright (c) 2020, A. Ovchinnikov, A. Pillay, G. Pogudin, T. Scanlon

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function _linear_compartment_model(graph, sinks)
    """
    Input: 
        - graph - graph of the model network represented via adjacency lists
        - sinks - the indices of nodes having a sink
    Output: the corresponding ODE object where each parameter a_ij is replaced
            with a_ij + b_ij * x_0, where x_0 is a constant input encoded as a constant output
    """
    n = length(graph)
    x_vars_names = ["x$i" for i in 0:n]
    edges_vars_names = Array{String, 1}()
    for i in 1:n
        for j in graph[i]
            push!(edges_vars_names, "a_$(i)_$(j)")
            push!(edges_vars_names, "b_$(i)_$(j)")
        end
    end
    for s in sinks
        push!(edges_vars_names, "a_$(s)_0")
        push!(edges_vars_names, "b_$(s)_0")
    end
    R, vars =
        Nemo.polynomial_ring(Nemo.QQ, vcat(x_vars_names, edges_vars_names, ["y1", "y2"]))
    x_vars = vars[2:(n + 1)]
    x0 = vars[1]
    equations =
        Dict{QQMPolyRingElem, Union{QQMPolyRingElem, Generic.Frac{QQMPolyRingElem}}}(
            x => R(0) for x in x_vars
        )
    equations[x0] = R(0)
    for i in 1:n
        for j in graph[i]
            rate = str_to_var("a_$(i)_$(j)", R) + str_to_var("b_$(i)_$(j)", R) * x0

            if i != j
                equations[x_vars[j]] += x_vars[i] * rate
                equations[x_vars[i]] -= x_vars[i] * rate
            else
                equations[x_vars[i]] -= x_vars[i] * rate
            end
        end
        if i in sinks
            rate = str_to_var("a_$(i)_0", R) + str_to_var("b_$(i)_0", R) * x0
            equations[x_vars[i]] += -x_vars[i] * rate
        end
    end
    return ODE{QQMPolyRingElem}(
        equations,
        Dict(vars[end] => x_vars[1], vars[end - 1] => x0),
        Array{QQMPolyRingElem, 1}(),
    )
end

#------------------------------------------------------------------------------

function bicycle(n)
    """
    Generates a bidirected cycle of length n
    """
    graph = []
    for i in 1:n
        prev = (i == 1) ? n : (i - 1)
        next = (i == n) ? 1 : i + 1
        push!(graph, [prev, next])
    end
    return graph
end

#------------------------------------------------------------------------------

function cycle(n)
    """
    Single directed cycle
    """
    graph = [[(i == n) ? 1 : (i + 1)] for i in 1:n]
    return graph
end

#------------------------------------------------------------------------------

function catenary(n)
    """
    Bidirected chain from 1 to n
    """
    graph = [[] for i in 1:n]
    for i in 1:n
        if i != 1
            push!(graph[i], i - 1)
        end
        if i != n
            push!(graph[i], i + 1)
        end
    end
    return graph
end

#------------------------------------------------------------------------------

function mammilary(n)
    """
    Bidirected 'star' with center at 1 and rays to 2, ..., n
    """
    graph = []
    push!(graph, [i for i in 2:n])
    for i in 2:n
        push!(graph, [1])
    end
    return graph
end

#------------------------------------------------------------------------------

###############################################################################

@testset "Assessing local identifiability (multiexperiment)" begin

    # checking bounds

    test_cases = [
        Dict(:name => "Cyclic", :graph => cycle, :bound => [3]),
        Dict(:name => "Catenary", :graph => catenary, :bound => [4, 5]),
        Dict(:name => "Mammilary", :graph => mammilary, :bound => [4, 5]),
    ]

    n_min = 3
    n_max = 8
    for case in test_cases
        for n in n_min:n_max
            model = _linear_compartment_model(case[:graph](n), [1])
            println(case[:name] * ", n = $n")
            @time result = assess_local_identifiability(model, type = :ME)
            correct = undef
            if n - n_min + 1 > length(case[:bound])
                correct = case[:bound][end]
            else
                correct = case[:bound][n - n_min + 1]
            end
            @test correct == result[2]
        end
    end

    # checking bounds and results

    test_cases = []

    ode = @ODEmodel(
        x0'(t) = -(a01 + a21) * x0(t) + a12 * x1(t),
        x1'(t) = a21 * x0(t) - a12 * x1(t),
        y(t) = x0(t)
    )
    funcs_to_test =
        [a01, a21, a12, a01 * a12, a01 + a12 + a21, (a01 + a12 + a21) // (a01 * a12)]
    correct = [false, false, false, true, true, true]
    push!(
        test_cases,
        Dict(
            :ode => ode,
            :funcs => funcs_to_test,
            :correct => (OrderedDict(funcs_to_test .=> correct), 1),
        ),
    )

    #--------------------------------------------------------------------------

    # Example 7.7 from https://arxiv.org/pdf/2011.10868.pdf
    ode = @ODEmodel(x0'(t) = 0, x1'(t) = x0(t) * x1(t) + mu1 * x0(t) + mu2, y(t) = x1(t))
    funcs_to_test = [mu1, mu2]
    correct = [true, true]
    push!(
        test_cases,
        Dict(
            :ode => ode,
            :funcs => funcs_to_test,
            :correct => (OrderedDict(funcs_to_test .=> correct), 2),
        ),
    )

    #--------------------------------------------------------------------------

    # Example 7.8 from https://arxiv.org/pdf/2011.10868.pdf
    ode = @ODEmodel(
        S'(t) = -b * S(t) * I(t) / N,
        E'(t) = b * S(t) * I(t) / N - nu * E(t),
        I'(t) = nu * E(t) - a * I(t),
        c'(t) = 0,
        y1(t) = c(t) * I(t) + d * E(t),
        y2(t) = c(t),
        y3(t) = N
    )
    funcs_to_test = [b, nu, d, a]
    correct = OrderedDict([b => true, nu => true, d => true, a => true])
    push!(test_cases, Dict(:ode => ode, :funcs => funcs_to_test, :correct => (correct, 1)))

    #--------------------------------------------------------------------------

    # example with 0 replicas required
    ode = @ODEmodel(x'(t) = a * z(t), z'(t) = a * z(t)^2, y(t) = x(t))
    funcs_to_test = [a]
    correct = OrderedDict([a => false])
    push!(test_cases, Dict(:ode => ode, :funcs => funcs_to_test, :correct => (correct, 0)))

    #--------------------------------------------------------------------------

    for case in test_cases
        result = assess_local_identifiability(
            case[:ode],
            funcs_to_check = case[:funcs],
            prob_threshold = 0.932,
            type = :ME,
        )
        @test result == case[:correct]
    end
end
