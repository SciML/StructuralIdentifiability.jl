# ------------------------------------------------------------------------------

function generator_var_change(generator, var::P, numer::P, denom::P) where P <: MPolyElem
    return IterTools.imap(
        point -> begin
            result = copy(point)
            result[var] = eval_at_dict(numer, point) // eval_at_dict(denom, point)
            return result
    end, Iterators.filter(point -> eval_at_dict(denom, point) != 0, generator)
    )
end

# ------------------------------------------------------------------------------

function diff_poly(poly::P, derivation::Dict{P,P}) where P <: MPolyElem
    return sum(derivative(poly, x) * xd for (x, xd) in derivation)
end

# ------------------------------------------------------------------------------

function generate_io_equation_problem(ode::ODE{P}) where P <: MPolyElem{<: FieldElem}
    dim_x = length(ode.x_vars)

    # Creating a ring
    old_vars = map(string, gens(ode.poly_ring))
    var_names = vcat(
        old_vars,
        [var_to_str(x) * "_dot" for x in ode.x_vars],
        [var_to_str(y) * "_$i" for i in 0:dim_x for y in ode.y_vars],
        [var_to_str(u) * "_$i" for i in 0:dim_x for u in ode.u_vars],
        ["rand_proj_var"]
    )
    ring, ring_vars = Nemo.PolynomialRing(base_ring(ode.poly_ring), var_names)

    # Definiting a (partial) derivation on it
    derivation = Dict{P,P}()
    for x in ode.x_vars
        derivation[switch_ring(x, ring)] = str_to_var(var_to_str(x) * "_dot", ring)
    end
    for i in 0:(dim_x - 1)
        for y in ode.y_vars
            derivation[str_to_var(var_to_str(y) * "_$i", ring)] = str_to_var(var_to_str(y) * "_$(i + 1)", ring)
        end
    end
    for u in ode.u_vars
        for i in 0:(dim_x - 1)
            derivation[str_to_var(var_to_str(u) * "_$i", ring)] = str_to_var(var_to_str(u) * "_$(i + 1)", ring)
        end
    end

    # Generate equations
    old_us = [parent_ring_change(u, ring) for u in ode.u_vars]
    new_us = [str_to_var(var_to_str(u) * "_0", ring) for u in ode.u_vars]
    function to_new_ring(p)
        return evaluate(parent_ring_change(ode.poly_ring(p), ring), old_us, new_us)
    end
    x_equations = Dict{P,P}()
    for x in ode.x_vars
        x_lifted = parent_ring_change(x, ring)
        num, den = map(to_new_ring, unpack_fraction(ode.x_equations[x]))
        x_equations[x_lifted] = den * derivation[x_lifted] - num
    end
    y_equations = Dict{P,P}()
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
    find_ioequations(ode, [var_change_policy=:default])

Finds the input-output equations of an ODE system
Input:
- `ode` - the ODE system
- `var_change_policy` - whether to perform automatic variable change, can be one of `:default`, `:yes`, `:no`

Output:
- a dictionary from "leaders" to the corresponding input-output equations
"""
function find_ioequations(
        ode::ODE{P}; 
        var_change_policy=:default
    ) where P <: MPolyElem{<: FieldElem}
    # Setting the var_change policy
    if (var_change_policy == :yes) || (var_change_policy == :default && length(ode.y_vars) < 3)
        auto_var_change = true
    elseif (var_change_policy == :no) || (var_change_policy == :default && length(ode.y_vars) >= 3)
        auto_var_change = false
    else
        @error "Unknown var_change policy $var_change_policy"
        return
    end

    # Initialization
    ring, derivation, x_equations, y_equations, point_generator = generate_io_equation_problem(ode)
    y_orders = Dict(y => 0 for y in keys(y_equations))
 
    while true        
        var_degs = [(y, [degree(eq, x) for x in keys(x_equations) if degree(eq, x) > 0]) for (y, eq) in y_equations]
        filter!(d -> length(d[2]) > 0, var_degs)
        if isempty(var_degs)
            break
        end
        @debug "Current degrees of io-equations $var_degs"
        @debug "Orders: $y_orders"
        @debug "Sizes: $(Dict(y => length(eq) for (y, eq) in y_equations))"

        # choosing the output to prolong
        outputs_with_scores = [
            (
                min(d[2]...) * length(y_equations[d[1]]),
                min(d[2]...),
                -count(x -> x == min(d[2]...), d[2]) + length(y_equations[d[1]]) // 30,
                length(y_equations[d[1]]),
                d[1]
            ) for d in var_degs 
        ]
        @debug "Scores: $outputs_with_scores"
        y_prolong = sort(outputs_with_scores)[1][end]
        y_orders[y_prolong] += 1
        @debug "Prolonging output $y_prolong"
        flush(stdout)

        # Calculate the Lie derivative of the io_relation
        @debug "Prolonging"
        flush(stdout)
        next_y_equation = diff_poly(y_equations[y_prolong], derivation)
        for (x, eq) in x_equations
            @debug "Eliminating the derivative of $x"
            flush(stdout)
            next_y_equation = eliminate_var(eq, next_y_equation, derivation[x], point_generator)
        end
        for (y, eq) in y_equations
            if y != y_prolong
                @debug "Eliminating the leader of the equation for $y"
                flush(stdout)
                # an ugly way of gettin the leader, to replace
                next_y_equation = eliminate_var(
                    next_y_equation, eq, 
                    str_to_var(var_to_str(y)[1:end - 2] * "_$(y_orders[y])", ring), 
                    point_generator
                )
            end
        end
        
        # Choose variable to eliminate
        var_degs_next = [(degree(y_equations[y_prolong], x), degree(next_y_equation, x), x) for x in keys(x_equations) if degree(y_equations[y_prolong], x) > 0]
        our_choice = sort(var_degs_next)[1]
        var_elim_deg, var_elim = our_choice[1], our_choice[3]
        
        @debug "Elimination of $var_elim, $(length(x_equations)) left; $(Dates.now())"
        flush(stdout)
        
        # Possible variable change for Axy + Bx + p(y) (x = var_elim)
        if auto_var_change && (var_elim_deg == 1)
            Ay_plus_B = coeff(y_equations[y_prolong], [var_elim], [1])
            for x in setdiff(keys(x_equations), [var_elim])
                if degree(Ay_plus_B, x) == 1                      
                    A, B = divrem(Ay_plus_B, x)
                    A, B = simplify_frac(A, B)
                    if isempty(filter!(v -> (v in keys(x_equations)), vars(A))) && (B != 0) # && (length(coeffs(A))==1) 
                        # variable change x_i' --> x_i' - (B/A)', x_i --> x_i - B/A
                        @debug "\t Applying variable change: $(x) --> $(x) - ( $B )/( $A ); $(Dates.now())"
                        flush(stdout)
                        dB = diff_poly(B, derivation)
                        dA = diff_poly(A, derivation)
                        numer_d, denom_d = simplify_frac(A * dB - dA * B, A * A)
                        # !!!push x_i_dot first, otherwise there will be problem with plugin!!!
                        point_generator = generator_var_change(
                            point_generator, 
                            derivation[x], 
                            derivation[x] * denom_d + numer_d, denom_d
                        )
                        point_generator = generator_var_change(point_generator, x, A * x + B, A)

                        # Change current system
                        @debug "Change in the system"
                        flush(stdout)
                        x_equations[x] = make_substitution(x_equations[x], derivation[x], denom_d * derivation[x] - numer_d, denom_d)
                        for xx in keys(x_equations)
                            @debug "\t Change in the equation for $xx"
                            flush(stdout)
                            x_equations[xx] = make_substitution(x_equations[xx], x, A * x - B, A)
                        end
                        @debug "Change in the outputs"
                        flush(stdout)
                        for y in keys(y_equations)
                            @debug "\t Change in the output $y"
                            flush(stdout)
                            y_equations[y] = make_substitution(y_equations[y], x, A * x - B, A)
                        end
                        @debug "\t Change in the prolonged equation"
                        flush(stdout)
                        next_y_equation = make_substitution(next_y_equation, x, A * x - B, A)
                        # recalibrate system
                        @debug "Unmixing the derivatives"
                        flush(stdout)
                        for xx in setdiff(keys(x_equations), [x])
                            @debug "\t Unmixing $xx"
                            flush(stdout)
                            x_equations[x] = eliminate_var(x_equations[x], x_equations[xx], derivation[xx], point_generator)
                        end
                        @debug "Change of variables performed"
                        flush(stdout)
                        break
                    end
                end
            end  
        end
   
        # Eliminate var_elim from the system
        delete!(x_equations, var_elim)
        @debug "Elimination in states"
        flush(stdout)
        for (x, eq) in x_equations
            @debug "\t Elimination in the equation for $x"
            flush(stdout)
            x_equations[x] = eliminate_var(eq, y_equations[y_prolong], var_elim, point_generator)
        end
        @debug "Elimination in y_equations"
        flush(stdout)
        for y in keys(y_equations)
            if y != y_prolong
                @debug "Elimination in the output $y"
                flush(stdout)
                y_equations[y] = eliminate_var(y_equations[y], y_equations[y_prolong], var_elim, point_generator)
            end
        end
        @debug "\t Elimination in the prolonged equation"
        flush(stdout)
        y_equations[y_prolong] = eliminate_var(y_equations[y_prolong], next_y_equation, var_elim, point_generator)
    end

    io_equations = Dict(str_to_var(var_to_str(y)[1:end - 2] * "_$(y_orders[y])", ring) => p for (y, p) in y_equations)

    @debug "Check whether the original projections are enough"
    if length(io_equations) == 1 || check_primality(io_equations) 
        @debug "The projections generate an ideal with a single components of highest dimension, returning"
        return io_equations
    end

    sampling_range = 5
    while true
        @debug "There are several components of the highest dimension, trying to isolate one"
        extra_var = str_to_var("rand_proj_var", ring)
        extra_var_rhs = sum(rand(1:sampling_range) * v for v in keys(io_equations))
        @debug "Extra projections: $extra_var_rhs"
        point_generator = generator_var_change(point_generator, extra_var, extra_var_rhs, one(ring))
        extra_poly = extra_var - extra_var_rhs
        for (y, io_eq) in io_equations
            extra_poly = eliminate_var(extra_poly, io_eq, y, point_generator)
        end
        extra_poly = evaluate(extra_poly, [extra_var], [extra_var_rhs])
        if check_primality(io_equations, [extra_poly])
            @debug "Single component of highest dimension isolated, returning"
            io_equations[extra_var] = extra_poly
            break
        end
        sampling_range = 2 * sampling_range
    end
    return io_equations
end
