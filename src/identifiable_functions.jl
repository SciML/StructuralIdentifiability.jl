
"""
    extract_identifiable_functions(io_equations, parameters)

For the io_equation and the list of all parameter variables, returns a set of generators of a field of all functions of parameters

Note: an experimental functionality at the moment, may fail and may be inefficient
"""
function extract_identifiable_functions(
    io_equations::Array{P, 1},
    parameters::Array{P, 1},
    known_functions::Array{P, 1},
) where {P <: MPolyElem{fmpq}}
    @debug "Extracting coefficients"
    flush(stdout)
    nonparameters = filter(
        v -> !(var_to_str(v) in map(var_to_str, parameters)),
        gens(parent(io_equations[1])),
    )
    coeff_lists = Array{Array{P, 1}, 1}()
    for eq in io_equations
        push!(coeff_lists, collect(values(extract_coefficients(eq, nonparameters))))
    end
    for f in known_functions
        push!(coeff_lists, [one(parent(f)), f])
    end
    for p in coeff_lists
        @debug sort(map(total_degree, p))
    end

    @debug "Resulting Coefficient List: $coeff_lists"

    return coeff_lists
end

#------------------------------------------------------------------------------

"""
    extract_identifiable_functions_raw(io_equations, parameters)

For the io_equation and the list of all parameter variables, returns a set of *raw* *generators of a field of all functions of parameters
"""
function extract_identifiable_functions_raw(
    io_equations::Array{P, 1},
    parameters::Array{P, 1},
) where {P <: MPolyElem{fmpq}}
    @debug "Extracting coefficients"
    flush(stdout)
    nonparameters = filter(
        v -> !(var_to_str(v) in map(var_to_str, parameters)),
        gens(parent(io_equations[1])),
    )
    result = []
    for eq in io_equations
        coeffs = sort(
            collect(values(extract_coefficients(eq, nonparameters))),
            by = total_degree,
        )
        append!(result, [c // first(coeffs) for c in coeffs[2:end]])
    end

    return result
end

#------------------------------------------------------------------------------

"""
    find_identifiable_functions(ode::ODE{<: MPolyElem{fmpq}}, p::Float64=0.99)

Input:
- `ode` - `ODE`-system
- `p` - probability of correctness

Output:
- returns a set of generators of the field of all functions of parameters

Find identifiable functions of parameters for a given `ode`.
"""
function find_identifiable_functions(
    ode::ODE{<:MPolyElem{fmpq}},
    p::Float64 = 0.99;
    simplify = true,
)
    @debug "Computing IO-equations"
    io_equations = find_ioequations(ode)
    global_result = check_identifiability(
        collect(values(io_equations)),
        ode.parameters,
        empty(ode.parameters),
        p,
    )
    known_params = Array{fmpq_mpoly, 1}()
    for (glob, p) in zip(global_result, ode.parameters)
        if glob
            push!(known_params, p)
        end
    end
    id_funcs = extract_identifiable_functions(
        collect(values(io_equations)),
        ode.parameters,
        known_params,
    )

    if simplify
        dideal = generate_diff_ideal(id_funcs)
        gb = ParamPunPam.paramgb(dideal)
        id_funcs = simplify_identifiable_functions(id_funcs, gb)
    end

    return id_funcs
end
