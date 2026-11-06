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

function suggest_variable_names(ode, new_states, new_params)
    @assert allunique(map(string, gens(parent(ode))))
    names = String[]

    for (i, def) in enumerate(new_states)
        num, den = unpack_fraction(def)
        name = nothing
        
        # x --> x
        if def in ode.x_vars
            name == nothing && (name = var_to_str(num))
        end

        # a*x+b --> x̃
        tildify(c::Char) = string(c, '\u0303')
        tildify(s) = tildify(first(s)) * s[nextind(s, firstindex(s)):end]
        states = intersect(vars(def), ode.x_vars)
        if length(states) == 1
            name_orig = var_to_str(states[1])
            name == nothing && (name = tildify(name_orig))
            # resolve collision a*x -> x̃ and b*x -> x̃
            j = 1
            while name in names
                name = chopsuffix(chopsuffix(name, "(t)"), string(j-1))
                name = string(name, j, "(t)")
                j += 1
            end
        end

        # f(x1,...,xn) --> X_i
        name == nothing && (name = "X$i(t)")

        push!(names, name)
    end

    for (i, def) in enumerate(new_params)
        num, den = unpack_fraction(def)
        name = nothing

        # x --> x
        if def in ode.parameters
            name == nothing && (name = var_to_str(num))
        end

        # f(x) --> ai
        name == nothing && (name = "a$i")

        push!(names, name)
    end

    append!(names, string.(ode.u_vars), string.(ode.y_vars))

    names
end

"""
    reparametrize_with_respect_to(ode, new_vars)

Reparametrizes the `ode` using the given new variables.

## Input

- `ode`: an ODE model.
- `new_vars`: a vector of new variables as fractions in `parent(ode)`.
"""
function reparametrize_with_respect_to(ode, new_vars)
    contains_states(poly::MPolyRingElem) = any(x -> degree(poly, x) > 0, ode.x_vars)
    contains_states(func) = contains_states(numerator(func)) || contains_states(denominator(func))
    new_vars_contains_states = filter(contains_states, new_vars)
    @info "Original variables:  $(length(ode.x_vars)) states and $(length(ode.parameters)) parameters"
    @debug """Original variables
    States: $(join(map(string, ode.x_vars), ", "))
    Params: $(join(map(string, ode.parameters), ", "))"""
    new_states = new_vars_contains_states
    new_params = setdiff(new_vars, new_vars_contains_states)
    @info "New macro variables: $(length(new_states)) states and $(length(new_params)) parameters"
    @debug """
    Reparametrizing with respect to:
    New states: $(join(map(string, new_states), ", "))
    New params: $(join(map(string, new_params), ", "))"""
    return reparametrize_with_respect_to(ode, new_states, new_params)
end

"""
    reparametrize_with_respect_to(ode, new_states, new_params)

Same as `reparametrize_with_respect_to(ode, new_vars)`,
but with the new states and new parameters specified separately.
"""
function reparametrize_with_respect_to(ode, new_states, new_params)
    @assert length(new_states) > 0
    poly_ring = base_ring(parent(first(new_states)))
    # Compute the new dynamics in terms of the original variables.
    # Paying attenton to the order..
    new_vector_field = vector_field_along(ode.x_equations, new_states)
    @debug "New vector field:\n$new_vector_field"
    new_states = collect(keys(new_vector_field))
    # Order by the first original state occurring in each new state.
    # TODO: might be too invasive if we do this here instead of the outer function
    sort!(new_states, by = f -> something(findfirst(in(vars(f)), ode.x_vars), length(ode.x_vars) + 1))
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
    @info "Constructing equations of new ODE"
    mem_time = @elapsed membership, new_dynamics_all, implicit_relations, new_vars =
        check_constructive_field_membership(
        RationalFunctionField(generating_funcs),
        to_be_reduced_funcs;
        tag_names = tag_names,
    )
    @info "Constructed equations of new ODE in $mem_time seconds"
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
    @info "Converting variable names to human-readable ones"
    new_names = suggest_variable_names(ode, new_states, new_params)
    @assert allunique(new_names)
    ring_output, _ = polynomial_ring(
        base_ring(ring_of_tags),
        new_names,
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
    new_vars = OrderedDict(
        parent_ring_change(var_old, ring_output, matching = :byindex) => new_vars[var_old] for
            var_old in gens(ring_of_tags)
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

function show_reparametrization(out, reparam)
    # eval(ode_repr(ode)) ≈ ode
    ode_repr(ode) = begin
        s = sprint(show, reparam.new_ode)
        s = join(map(line -> string(" "^2, line), split(rstrip(s, '\n'), "\n")), ",\n")
        string("ode = @ODEmodel(\n", s, "\n)")
    end

    println(out, "\n", "• Reparametrized ODE:\n")
    println(out, ode_repr(reparam.new_ode), "\n")

    println(out, "• Definitions of new variables:", "\n")
    println(out, join(map(p -> replace(string(" "^2, p), "=>" => ":="), collect(reparam.new_vars)), "\n"), "\n")

    if !isempty(reparam.implicit_relations)
        participants = mapreduce(vars, union, reparam.implicit_relations)
        println(out, "• These variables are not independent: $(join(map(string, participants), ", "))", "\n")
    end

    println(out, "• For more info, run\n\njulia> @doc reparametrize_global\n")
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
- `cmp`: A comparator for rational functions. For two rational functions `f1`
    and `f2`, `cmp(f1, f2)` should be `true` if `f1` is simpler than `f2`.
    By default, `StructuralIdentifiability.cmp_default(ode)` is used.

## Example

```jldoctest
julia> using StructuralIdentifiability

julia> using Logging

julia> Logging.disable_logging(Logging.Info);

julia> ode = @ODEmodel(
           x1'(t) = a * x1(t) - b*x1(t)*x2(t),
           x2'(t) = -c * x2(t) + d*x1(t)*x2(t),
           y(t) = x1(t)
       );

julia> Logging.disable_logging(Logging.BelowMinLevel);

julia> new_ode, new_vars, relations = reparametrize_global(
           ode; loglevel = Logging.Error
       );
```

Then, we have the following:

```
# new_ode
x1'(t) = a*x1(t) - x1(t)*x̃2(t)
x̃2'(t) = d*x1(t)*x̃2(t) - c*x̃2(t)
y(t) = x1(t)

# new_vars
Dict{Nemo.QQMPolyRingElem, AbstractAlgebra.Generic.FracFieldElem{Nemo.QQMPolyRingElem}} with 6 entries:
  x1(t) => x1(t)
  x̃2(t) => b*x2(t)
  y(t)  => y(t)
  d     => d
  a     => a
  c     => c
```

Notice that the `new_ode` is fully identifiabile, and has `1` less parameter
compared to the original one.
"""
function reparametrize_global(
        ode::ODE{P};
        prob_threshold = 0.99,
        seed = 42,
        cmp = cmp_default(ode),
        loglevel = Logging.Info,
    ) where {P}
    restart_logging(loglevel = loglevel)
    return with_logger(_si_logger[]) do
        return _reparametrize_global(ode, prob_threshold = prob_threshold, seed = seed, cmp = cmp)
    end
end

function _reparametrize_global(ode::ODE{P}; prob_threshold = 0.99, seed = 42, cmp = cmp_default(ode)) where {P}
    Random.seed!(seed)
    id_funcs = _find_identifiable_functions(
        ode,
        with_states = true,
        simplify = :strong,
        prob_threshold = prob_threshold,
        cmp = cmp,
    )
    ode_ring = parent(ode)
    @assert base_ring(parent(first(id_funcs))) == ode_ring
    @info "Constructing a new parametrization"
    new_vector_field, new_inputs, new_outputs, new_vars, implicit_relations =
        reparametrize_with_respect_to(ode, id_funcs)
    new_ring = parent(first(keys(new_vector_field)))
    # Fix the order of variables in the new ODE
    new_x_vars = filter(v -> haskey(new_vector_field, v), gens(new_ring))
    new_y_vars = filter(v -> haskey(new_outputs, v), gens(new_ring))
    new_ode = ODE{P}(new_x_vars, new_y_vars, new_vector_field, new_outputs, new_inputs)
    reparam = (new_ode = new_ode, new_vars = new_vars, implicit_relations = implicit_relations)
    show_reparametrization(stdout, reparam)
    return reparam
end
