"""
    vector_field_along(derivation, directions)

Returns the vector field obtained by applying `derivation` to each element of
`directions`.
"""
function vector_field_along(derivation::Dict{T, U}, directions::AbstractVector) where {T, U}
    new_vector_field = Dict{
        AbstractAlgebra.Generic.FracFieldElem{T},
        AbstractAlgebra.Generic.FracFieldElem{T},
    }()
    for func in directions
        df = diff_frac(func, derivation)
        new_vector_field[func] = df
    end
    return new_vector_field
end

function default_variable_names(new_states, new_params)
    (
        states = map(i -> "X$i", 1:length(new_states)),
        params = map(i -> "a$i", 1:length(new_params)),
    )
end

"""
    reparametrize_with_respect_to(ode, new_states, new_params)

Reparametrizes the `ode` using the given fractional states and parameters.

## Input

- `ode`: an ODE model.
- `new_states`: a vector of new states as fractions in `parent(ode)`.
- `new_params`: a vector of new parameters as fractions in `parent(ode)`.
"""
function reparametrize_with_respect_to(
    ode::ODE{P},
    new_states,
    new_params;
    new_variable_names = default_variable_names(new_states, new_params),
) where {P}
    @assert length(new_states) > 0
    poly_ring = base_ring(parent(first(new_states)))
    # Compute the new dynamics in terms of the original variables.
    # Paying attenton to the order..
    new_vector_field = vector_field_along(ode.x_equations, new_states)
    @debug "New vector field:\n$new_vector_field"
    new_states = collect(keys(new_vector_field))
    new_dynamics = [new_vector_field[new_state] for new_state in new_states]
    # Express the new dynamics in terms of new states and new parameters.
    outputs = [ode.y_equations[output] for output in ode.y_vars]
    generating_funcs = vcat(
        new_states,
        new_params,
        ode.u_vars .// one(poly_ring),
        ode.y_vars .// one(poly_ring),
    )
    to_be_reduced_funcs = vcat(new_dynamics, outputs .// one(poly_ring))
    n_active_generators =
        (length(generating_funcs) - length(ode.u_vars) - length(ode.y_vars))
    tag_names = vcat(
        gen_tag_names(n_active_generators, "Internal"),
        gen_tag_names(length(ode.u_vars), "Input"),
        gen_tag_names(length(ode.y_vars), "Output"),
    )
    @debug """
    Tag names: 
    $tag_names
    Generating functions:
    $generating_funcs
    Functions to be reduced:
    $to_be_reduced_funcs
    """
    membership, new_dynamics_all, implicit_relations, new_vars =
        check_constructive_field_membership(
            RationalFunctionField(generating_funcs),
            to_be_reduced_funcs;
            tag_names = tag_names,
        )
    @assert all(membership)
    ring_of_tags = parent(first(keys(new_vars)))
    tags = gens(ring_of_tags)
    tag_inputs = tags[(n_active_generators + 1):(end - length(ode.y_vars))]
    tag_outputs = tags[(end - length(ode.y_vars) + 1):end]
    new_dynamics_states = new_dynamics_all[1:length(new_states)]
    new_dynamics_outputs = new_dynamics_all[(length(new_states) + 1):end]
    new_outputs = Dict(
        output => dynamic for (output, dynamic) in zip(tag_outputs, new_dynamics_outputs)
    )
    # Old inputs map one to one to new inputs.
    new_inputs = empty(tags)
    if !isempty(ode.u_vars)
        new_inputs = tag_inputs
    end
    @debug """
    New state dynamics:
    $new_dynamics_states
    New output dynamics:
    $new_dynamics_outputs
    New inputs:
    $new_inputs"""
    # Construct the new vector field.
    new_vars_vector_field = empty(ode.x_equations)
    for i in 1:length(new_states)
        state = tags[i]
        new_vars_vector_field[state] = new_dynamics_states[i]
    end
    @debug "Converting variable names to human-readable ones"
    @assert length(new_variable_names.states) == length(new_states)
    @assert length(new_variable_names.params) == length(new_params)
    input_variable_names = map(i -> "u$i(t)", 1:length(tag_inputs))
    output_variable_names = map(i -> "y$i(t)", 1:length(tag_outputs))
    all_variable_names = vcat(
        new_variable_names.states,
        new_variable_names.params,
        input_variable_names,
        output_variable_names,
    )
    ring_output, _ = polynomial_ring(
        base_ring(ring_of_tags),
        all_variable_names,
        internal_ordering = Nemo.internal_ordering(ring_of_tags),
    )
    new_vars_vector_field = Dict(
        parent_ring_change(var_old, ring_output, matching = :byindex) =>
            parent_ring_change(var_expr, ring_output, matching = :byindex) for
        (var_old, var_expr) in new_vars_vector_field
    )
    new_inputs = map(
        var_old -> parent_ring_change(var_old, ring_output, matching = :byindex),
        new_inputs,
    )
    new_outputs = Dict(
        parent_ring_change(var_old, ring_output, matching = :byindex) =>
            parent_ring_change(var_expr, ring_output, matching = :byindex) for
        (var_old, var_expr) in new_outputs
    )
    new_vars = Dict(
        parent_ring_change(var_old, ring_output, matching = :byindex) => var_expr for
        (var_old, var_expr) in new_vars
    )
    implicit_relations = map(
        var_old -> parent_ring_change(var_old, ring_output, matching = :byindex),
        implicit_relations,
    )
    @assert parent(first(keys(new_vars_vector_field))) ==
            base_ring(parent(first(values(new_vars_vector_field)))) ==
            parent(first(keys(new_outputs))) ==
            base_ring(parent(first(values(new_outputs)))) ==
            parent(first(keys(new_vars)))
    @assert base_ring(parent(first(values(new_vars)))) == parent(ode)
    new_vars_vector_field, new_inputs, new_outputs, new_vars, implicit_relations
end

"""
    reparametrize_global(ode, options...)

Finds a reparametrization of `ode` in terms of globally identifiabile functions.

Returns a tuple (`new_ode`, `new_vars`, `implicit_relations`), such that:
- `new_ode` is the reparametrized ODE system
- `new_vars` is a mapping from the new variables to the original ones
- `relations` is the array of implicit algebraic relations among `new_vars`.
  Usually, `relations` is empty.

## Options

The function accepts the following optional arguments.

- `seed`: A float in the range from 0 to 1, random seed (default is `seed = 42`). 
- `prob_threshold`: The probability of correctness (default is `prob_threshold = 0.99`).

## Example

```jldoctest
using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = a * x1(t) - b*x1(t)*x2(t),
    x2'(t) = -c * x2(t) + d*x1(t)*x2(t),
    y(t) = x1(t)
)

new_ode, new_vars, relations = reparametrize_global(ode)
```

Then, we have the following:

```
# new_ode
X2'(t) = X1(t)*X2(t)*a2 - X2(t)*a1
X1'(t) = -X1(t)*X2(t) + X1(t)*a3
y1(t) = X1(t)

# new_vars
Dict{Nemo.QQMPolyRingElem, AbstractAlgebra.Generic.FracFieldElem{Nemo.QQMPolyRingElem}} with 6 entries:
  X2 => b*x2
  y1 => y
  X1 => x1
  a2 => d
  a3 => a
  a1 => c
```

Notice that the `new_ode` is fully identifiabile, and has `1` less parameter
compared to the original one.
"""
function reparametrize_global(
    ode::ODE{P};
    prob_threshold = 0.99,
    seed = 42,
    loglevel = Logging.Info,
) where {P}
    restart_logging(loglevel = loglevel)
    with_logger(_si_logger[]) do
        return _reparametrize_global(ode, prob_threshold = prob_threshold, seed = seed)
    end
end

function _reparametrize_global(ode::ODE{P}; prob_threshold = 0.99, seed = 42) where {P}
    Random.seed!(seed)
    id_funcs = find_identifiable_functions(
        ode,
        with_states = true,
        simplify = :strong,
        prob_threshold = prob_threshold,
    )
    ode_ring = parent(ode)
    @assert base_ring(parent(first(id_funcs))) == ode_ring
    @info "Constructing a new parametrization"
    contains_states(poly::MPolyRingElem) = any(x -> degree(poly, x) > 0, ode.x_vars)
    contains_states(func) =
        contains_states(numerator(func)) || contains_states(denominator(func))
    id_funcs_contains_states = filter(contains_states, id_funcs)
    @info """
    Original states: $(ode.x_vars)
    Original params: $(ode.parameters)
    Identifiable and contain states: $(id_funcs_contains_states)"""
    new_states = id_funcs_contains_states
    new_params = setdiff(id_funcs, id_funcs_contains_states)
    @info """
    Reparametrizing with respect to:
    New states: $new_states
    New params: $new_params"""
    new_vector_field, new_inputs, new_outputs, new_vars, implicit_relations =
        reparametrize_with_respect_to(ode, new_states, new_params)
    new_ring = parent(first(keys(new_vector_field)))
    new_ode = ODE{P}(new_vector_field, new_outputs, new_inputs)
    return (new_ode = new_ode, new_vars = new_vars, implicit_relations = implicit_relations)
end

function reparametrize_interactive(
    ode::ODE{P};
    prob_threshold = 0.99,
    seed = 42,
    loglevel = Logging.Info,
    input::IO = stdin,
    output::IO = stdout,
) where {P}
    restart_logging(loglevel = loglevel, stream = output)
    with_logger(_si_logger[]) do
        return _reparametrize_interactive(ode, prob_threshold, seed, input, output)
    end
end

function _reparametrize_interactive(
    ode::ODE{P},
    prob_threshold,
    seed,
    input,
    output,
) where {P}
    Random.seed!(seed)
    id_funcs = find_identifiable_functions(
        ode,
        with_states = true,
        simplify = :strong,
        prob_threshold = prob_threshold,
        return_all = true,
    )
    id_funcs_simple = find_identifiable_functions(
        ode,
        with_states = true,
        simplify = :strong,
        prob_threshold = prob_threshold,
    )
    state = Dict(
        :counter => 0,
        :id_funcs => id_funcs,
        :chosen_funcs => (states = empty(id_funcs), params = empty(id_funcs)),
        :variable_names => (states = Vector{String}(), params = Vector{String}()),
    )
    contains_states(poly::MPolyRingElem) = any(x -> degree(poly, x) > 0, ode.x_vars)
    contains_states(func) =
        contains_states(numerator(func)) || contains_states(denominator(func))
    function print_header(state)
        println(output, "\n$(state[:counter]). Info.")
        println(output, "  Original states: $(join(string.(ode.x_vars), ", "))")
        println(output, "  Original parameters: $(join(string.(ode.parameters), ", "))")
        println(output, "  Identifiable functions: $(join(string.(id_funcs_simple), ", "))")
    end
    function print_state(state)
        counter, chosen_funcs, variable_names =
            state[:counter], state[:chosen_funcs], state[:variable_names]
        if isempty(chosen_funcs.states) && isempty(chosen_funcs.params)
            return
        end
        println(output, "\n$counter. Current selection:")
        for (names, funcs) in [
            (variable_names.states, chosen_funcs.states),
            (variable_names.params, chosen_funcs.params),
        ]
            for (name, func) in zip(names, funcs)
                println(output, "  ", name, " := ", func)
            end
        end
    end
    function make_choice(state)
        counter, id_funcs = state[:counter], state[:id_funcs]
        terminal =
            REPL.TerminalMenus.default_terminal(in = input, out = output, err = output)
        menu = MultiSelectMenu(vcat("Enter a custom function", map(string, id_funcs)))
        choice = request(
            terminal,
            "\n$counter. Select identifiable function(s) for reparametrization.",
            menu,
        )
        if 1 in choice # a custom function
            funcs = empty(id_funcs)
            while true
                varnames = map(
                    f -> chopsuffix(f, "(t)"),
                    string.(vcat(ode.x_vars, ode.parameters)),
                )
                res = Base.prompt(
                    input,
                    output,
                    "\n$counter. Enter a rational function in the variables: $(join(varnames, ", "))\n",
                )
                func = nothing
                try
                    func = myeval(
                        Meta.parse(res),
                        Dict(Symbol.(varnames) .=> vcat(ode.x_vars, ode.parameters)),
                    )
                catch e
                    @info "" e
                    printstyled(
                        output,
                        "\n  ==> Error when parsing $res. Trying again..\n",
                        bold = true,
                    )
                    continue
                end
                ffring = fraction_field(parent(ode))
                funcs = [ffring(func)]
                if all(
                    field_contains(
                        RationalFunctionField(id_funcs_simple),
                        funcs,
                        prob_threshold,
                    ),
                )
                    break
                else
                    printstyled(
                        output,
                        "\n  ==> The given function $(funcs[1]) is not identifiable. Trying again..\n",
                        bold = true,
                    )
                    continue
                end
            end
        else
            funcs = id_funcs[sort(collect(choice)) .- 1]
        end
        funcs
    end
    function query_names(state, funcs)
        counter, chosen_funcs, variable_names =
            state[:counter], state[:chosen_funcs], state[:variable_names]
        new_states = filter(contains_states, funcs)
        new_params = setdiff(funcs, new_states)
        idx_states, idx_params = length(chosen_funcs.states), length(chosen_funcs.params)
        append!(chosen_funcs.states, new_states)
        append!(chosen_funcs.params, new_params)
        default_names = default_variable_names(chosen_funcs.states, chosen_funcs.params)
        for (kind, vars, defaults, new_vars, idx) in [
            ("state", new_states, default_names.states, variable_names.states, idx_states),
            (
                "parameter",
                new_params,
                default_names.params,
                variable_names.params,
                idx_params,
            ),
        ]
            for (i, var) in enumerate(vars)
                default = defaults[idx + i]
                res = Base.prompt(
                    input,
                    output,
                    "\n$counter. Enter a name for the new $kind: $var. Leave empty for default: $default.\n",
                )
                if isnothing(res) || (res isa String && isempty(strip(res)))
                    res = default
                end
                push!(new_vars, res)
                printstyled(output, "  ==> New variable: $res := $var\n", bold = true)
            end
        end
    end
    function try_to_reparametrize(state)
        counter, chosen_funcs, variable_names, id_funcs =
            state[:counter], state[:chosen_funcs], state[:variable_names], state[:id_funcs]
        rff = RationalFunctionField(vcat(chosen_funcs.states, chosen_funcs.params))
        state[:id_funcs] = id_funcs[.! field_contains(rff, id_funcs, prob_threshold)]
        if isempty(chosen_funcs.states)
            printstyled(
                output,
                "\n  ==> Please select at least one new state in order to reparametrize.\n",
                bold = true,
            )
            return nothing
        end
        try
            new_vector_field, new_inputs, new_outputs, new_vars, implicit_relations =
                reparametrize_with_respect_to(
                    ode,
                    chosen_funcs.states,
                    chosen_funcs.params,
                    new_variable_names = variable_names,
                )
            new_ode = ODE{P}(new_vector_field, new_outputs, new_inputs)
            return (
                new_ode = new_ode,
                new_vars = new_vars,
                implicit_relations = implicit_relations,
            )
        catch e
            printstyled(
                output,
                "\n  ==> Selected functions is not enough to reparametrize. Please select more.\n",
                bold = true,
            )
        end
        return nothing
    end

    print_header(state)
    while true
        state[:counter] += 1
        print_state(state)
        funcs = make_choice(state)
        query_names(state, funcs)
        result = try_to_reparametrize(state)
        !(result == nothing) && return result
    end
end
