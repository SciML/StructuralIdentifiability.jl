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

"""
    reparametrize_with_respect_to(ode, new_states, new_params)

Reparametrizes the `ode` using the given fractional states and parameters.

## Input

- `ode`: an ODE model.
- `new_states`: a vector of new states as fractions in `parent(ode)`.
- `new_params`: a vector of new parameters as fractions in `parent(ode)`.
"""
function reparametrize_with_respect_to(ode, new_states, new_params)
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
    @info """
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
    @info """
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
    @info "Converting variable names to human-readable ones"
    internal_variable_names = map(i -> "X$i(t)", 1:length(new_states))
    parameter_variable_names = map(i -> "a$i", 1:length(new_params))
    input_variable_names = map(i -> "u$i(t)", 1:length(tag_inputs))
    output_variable_names = map(i -> "y$i(t)", 1:length(tag_outputs))
    all_variable_names = vcat(
        internal_variable_names,
        parameter_variable_names,
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
    return new_vars_vector_field, new_inputs, new_outputs, new_vars, implicit_relations
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
- `loglevel`: the level of logs to be displayed (default is `Logging.Info`).

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
    return with_logger(_si_logger[]) do
        return _reparametrize_global(ode, prob_threshold = prob_threshold, seed = seed)
    end
end

function _reparametrize_global(ode::ODE{P}; prob_threshold = 0.99, seed = 42) where {P}
    Random.seed!(seed)
    id_funcs = _find_identifiable_functions(
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
