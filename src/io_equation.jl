function generator_var_change(
        generator,
        var::P,
        numer::P,
        denom::P,
    ) where {P <: MPolyRingElem}
    return IterTools.imap(
        point -> begin
            result = copy(point)
            result[var] = eval_at_dict(numer, point) // eval_at_dict(denom, point)
            return result
        end,
        Iterators.filter(point -> eval_at_dict(denom, point) != 0, generator),
    )
end

# ------------------------------------------------------------------------------

function diff_poly(poly::P, derivation::Dict{P, T}) where {P <: MPolyRingElem, T}
    return sum(derivative(poly, x) * xd for (x, xd) in derivation)
end

function diff_frac(frac::Generic.FracFieldElem{F}, derivation::Dict{P, T}) where {F, P, T}
    num, den = unpack_fraction(frac)
    numd = den * diff_poly(num, derivation) - num * diff_poly(den, derivation)
    dend = den^2
    return numd // dend
end

# ------------------------------------------------------------------------------

function generate_io_equation_problem(ode::ODE{P}) where {P <: MPolyRingElem{<:FieldElem}}
    dim_x = length(ode.x_vars)

    # Creating a ring
    old_vars = map(string, gens(ode.poly_ring))
    var_names = vcat(
        old_vars,
        [var_to_str(x) * "_dot" for x in ode.x_vars],
        [var_to_str(y) * "_$i" for i in 0:dim_x for y in ode.y_vars],
        [var_to_str(u) * "_$i" for i in 0:dim_x for u in ode.u_vars],
    )
    ring, ring_vars = Nemo.polynomial_ring(base_ring(ode.poly_ring), var_names)

    # Defining a (partial) derivation on it
    derivation = Dict{P, P}()
    for x in ode.x_vars
        derivation[switch_ring(x, ring)] = str_to_var(var_to_str(x) * "_dot", ring)
    end
    for i in 0:(dim_x - 1)
        for y in ode.y_vars
            derivation[str_to_var(var_to_str(y) * "_$i", ring)] =
                str_to_var(var_to_str(y) * "_$(i + 1)", ring)
        end
    end
    for u in ode.u_vars
        for i in 0:(dim_x - 1)
            derivation[str_to_var(var_to_str(u) * "_$i", ring)] =
                str_to_var(var_to_str(u) * "_$(i + 1)", ring)
        end
    end

    # Generate equations
    old_us = [parent_ring_change(u, ring) for u in ode.u_vars]
    new_us = [str_to_var(var_to_str(u) * "_0", ring) for u in ode.u_vars]
    function to_new_ring(p)
        return evaluate(parent_ring_change(ode.poly_ring(p), ring), old_us, new_us)
    end
    x_equations = Dict{P, P}()
    for x in ode.x_vars
        x_lifted = parent_ring_change(x, ring)
        num, den = map(to_new_ring, unpack_fraction(ode.x_equations[x]))
        x_equations[x_lifted] = den * derivation[x_lifted] - num
    end
    y_equations = Dict{P, P}()
    for y in ode.y_vars
        y_lifted = str_to_var(var_to_str(y) * "_0", ring)
        g_num, g_den = map(to_new_ring, unpack_fraction(ode.y_equations[y]))
        y_equations[y_lifted] = g_den * y_lifted - g_num
    end

    generic_point_generator = ODEPointGenerator{P}(ode, ring)

    return (ring, derivation, x_equations, y_equations, generic_point_generator)
end

# ------------------------------------------------------------------------------

"""
    find_ioprojections(ode, auto_var_change)

Finds the input-output projections of an ODE system
Input:
- `ode` - the ODE system
- `auto_var_change` - whether to perform automatic variable change

Output:
- a dictionary from “leaders” to the corresponding input-output equations
- generic point generator for the model (including the derivatives; mostly for testing)
"""
@timeit _to function find_ioprojections(
        ode::ODE{P},
        auto_var_change::Bool,
    ) where {P <: MPolyRingElem{<:FieldElem}}
    # Initialization
    ring, derivation, x_equations, y_equations, point_generator =
        generate_io_equation_problem(ode)
    y_orders = Dict(y => 0 for y in keys(y_equations))

    while true
        var_degs = [
            (y, [Nemo.degree(eq, x) for x in keys(x_equations) if Nemo.degree(eq, x) > 0]) for (y, eq) in y_equations
        ]
        filter!(d -> length(d[2]) > 0, var_degs)
        if isempty(var_degs)
            break
        end
        @debug "Current degrees of io-equations $var_degs"
        @debug "Orders: $y_orders"
        @debug "Sizes: $(Dict(y => length(eq) for (y, eq) in y_equations))"
        flush(_si_logger[].stream)

        # choosing the output to prolong
        outputs_with_scores = [
            (
                min(d[2]...) * length(y_equations[d[1]]),
                min(d[2]...),
                -count(x -> x == min(d[2]...), d[2]) + length(y_equations[d[1]]) // 30,
                length(y_equations[d[1]]),
                d[1],
            ) for d in var_degs
        ]
        @debug "Scores: $outputs_with_scores"
        y_prolong = sort(outputs_with_scores)[1][end]
        y_orders[y_prolong] += 1
        @debug "Prolonging output $y_prolong"

        # Calculate the Lie derivative of the io_relation
        @debug "Prolonging"
        next_y_equation = diff_poly(y_equations[y_prolong], derivation)
        for (x, eq) in x_equations
            @debug "Eliminating the derivative of $x"
            next_y_equation =
                eliminate_var(eq, next_y_equation, derivation[x], point_generator)
        end
        for (y, eq) in y_equations
            if y != y_prolong
                @debug "Eliminating the leader of the equation for $y"
                # an ugly way of getting the leader, to replace
                next_y_equation = eliminate_var(
                    next_y_equation,
                    eq,
                    str_to_var(var_to_str(y)[1:(end - 1)] * "$(y_orders[y])", ring),
                    point_generator,
                )
            end
        end

        # Choose variable to eliminate
        var_degs_next = [
            (Nemo.degree(y_equations[y_prolong], x), Nemo.degree(next_y_equation, x), x) for x in keys(x_equations) if Nemo.degree(y_equations[y_prolong], x) > 0
        ]
        @debug "Choosing the variable to eliminate, scores: $(sort(var_degs_next))"
        our_choice = sort(var_degs_next)[1]
        @debug "Chosen: $our_choice"
        var_elim_deg, var_elim = our_choice[1], our_choice[3]

        @debug "Elimination of $var_elim, $(length(x_equations)) left"

        # Possible variable change for Axy + Bx + p(y) (x = var_elim)
        if auto_var_change && (var_elim_deg == 1)
            Ay_plus_B = coeff(y_equations[y_prolong], [var_elim], [1])
            for x in setdiff(keys(x_equations), [var_elim])
                if Nemo.degree(Ay_plus_B, x) == 1
                    A, B = divrem(Ay_plus_B, x)
                    A, B = simplify_frac(A, B)
                    if isempty(filter!(v -> (v in keys(x_equations)), vars(A))) && (B != 0) # && (length(coeffs(A))==1)
                        # variable change x_i' --> x_i' - (B/A)', x_i --> x_i - B/A
                        @debug "\t Applying variable change: $(x) --> $(x) - ( $B )/( $A )"
                        dB = diff_poly(B, derivation)
                        dA = diff_poly(A, derivation)
                        numer_d, denom_d = simplify_frac(A * dB - dA * B, A * A)
                        # !!!push x_i_dot first, otherwise there will be problem with plugin!!!
                        point_generator = generator_var_change(
                            point_generator,
                            derivation[x],
                            derivation[x] * denom_d + numer_d,
                            denom_d,
                        )
                        point_generator =
                            generator_var_change(point_generator, x, A * x + B, A)

                        # Change current system
                        @debug "Change in the system"
                        x_equations[x] = make_substitution(
                            x_equations[x],
                            derivation[x],
                            denom_d * derivation[x] - numer_d,
                            denom_d,
                        )
                        for xx in keys(x_equations)
                            @debug "\t Change in the equation for $xx"
                            x_equations[xx] =
                                make_substitution(x_equations[xx], x, A * x - B, A)
                        end
                        @debug "Change in the outputs"
                        for y in keys(y_equations)
                            @debug "\t Change in the output $y"
                            y_equations[y] =
                                make_substitution(y_equations[y], x, A * x - B, A)
                        end
                        @debug "\t Change in the prolonged equation"
                        next_y_equation =
                            make_substitution(next_y_equation, x, A * x - B, A)
                        # recalibrate system
                        @debug "Unmixing the derivatives"
                        for xx in setdiff(keys(x_equations), [x])
                            @debug "\t Unmixing $xx"
                            x_equations[x] = eliminate_var(
                                x_equations[x],
                                x_equations[xx],
                                derivation[xx],
                                point_generator,
                            )
                        end
                        break
                    end
                end
            end
        end

        # Eliminate var_elim from the system
        delete!(x_equations, var_elim)
        @debug "Elimination in states"
        for (x, eq) in x_equations
            @debug "\t Elimination in the equation for $x"
            x_equations[x] =
                eliminate_var(eq, y_equations[y_prolong], var_elim, point_generator)
        end
        @debug "Elimination in y_equations"
        for y in keys(y_equations)
            if y != y_prolong
                @debug "Elimination in the output $y"
                y_equations[y] = eliminate_var(
                    y_equations[y],
                    y_equations[y_prolong],
                    var_elim,
                    point_generator,
                )
            end
        end
        @debug "\t Elimination in the prolonged equation"
        y_equations[y_prolong] = eliminate_var(
            y_equations[y_prolong],
            next_y_equation,
            var_elim,
            point_generator,
        )
    end

    io_projections = Dict(
        str_to_var(var_to_str(y)[1:(end - 1)] * "$(y_orders[y])", ring) => p for
            (y, p) in y_equations
    )

    return io_projections, point_generator
end

# ------------------------------------------------------------------------------

"""
    find_ioequations(ode, [var_change_policy=:default])

Finds the input-output equations of an ODE system
Input:
- `ode` - the ODE system
- `var_change_policy` - whether to perform automatic variable change, can be one of `:default`, `:yes`, `:no`
- `loglevel` - logging level (default: Logging.Info)

Output: a dictionary from “leaders” to the corresponding input-output equations
"""
@timeit _to function find_ioequations(
        ode::ODE{P};
        var_change_policy = :default,
        loglevel = Logging.Info,
    ) where {P <: MPolyRingElem{<:FieldElem}}
    restart_logging(loglevel = loglevel)
    reset_timings()
    with_logger(_si_logger[]) do
        ioeqs, is_prime = _find_ioequations(ode, var_change_policy = var_change_policy)
        if !is_prime
            @warn "The computed set of IO-projections may not define a prime ideal (in particular, the coefficients may not generate all identifiable functions)"
        end
        return ioeqs
    end
end

@timeit _to function _find_ioequations(
        ode::ODE{P};
        var_change_policy = :default,
    ) where {P <: MPolyRingElem{<:FieldElem}}
    # Setting the var_change policy
    if (var_change_policy == :yes) ||
            (var_change_policy == :default && length(ode.y_vars) < 3)
        auto_var_change = true
    elseif (var_change_policy == :no) ||
            (var_change_policy == :default && length(ode.y_vars) >= 3)
        auto_var_change = false
    else
        @error "Unknown var_change policy $var_change_policy"
        return
    end

    io_projections, _ = find_ioprojections(ode, auto_var_change)
    ring = parent(first(values(io_projections)))

    @debug "Sizes $(Dict(k => length(v) for (k, v) in io_projections))"
    @debug "Check whether the original projections are enough"
    if length(io_projections) == 1 || check_primality(io_projections)
        @debug "The projections generate an ideal with a single components of highest dimension"
        return io_projections, true
    end
    @debug "The projections may generate a non-prime ideal"
    return io_projections, false
end
