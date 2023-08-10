
function reparametrize(ode::ODE{T}) where {T <: MPolyElem{fmpq}}
    runtime_start = time_ns()
    @info "Finding identifiabile functions in states"
    id_funcs = find_identifiable_functions(ode, with_states = true, strategy = (:hybrid,))
    @info "Constructing a new parametrization"
    derivation = ode.x_equations
    contains_states(poly::MPolyElem) = any(x -> degree(poly, x) > 0, ode.x_vars)
    contains_states(func) =
        contains_states(numerator(func)) || contains_states(denominator(func))
    id_funcs_contains_states = filter(contains_states, id_funcs)
    @info """
    Original variables: $(ode.x_vars)
    Identifiable states: $(id_funcs_contains_states)"""
    new_vector_field = Dict()
    for id_func in id_funcs_contains_states
        df = diff_frac(id_func, derivation)
        new_vector_field[id_func] = df
    end
    new_vector_field
end