using Dates
using IterTools
using Logging
using Oscar

include("util.jl")
include("elimination.jl")
include("ODE.jl")

#------------------------------------------------------------------------------

function generator_var_change(generator, var, numer, denom)
    return IterTools.imap(
        point -> begin
            result = copy(point)
            result[var] = eval_at_dict(numer, point) // eval_at_dict(denom, point)
            return result
        end,
        Iterators.filter(point -> eval_at_dict(denom, point) != 0, generator)
    )
end

#------------------------------------------------------------------------------

function diff_poly(poly, derivation)
    return sum([derivative(poly, x) * xd for (x, xd) in derivation])
end

#------------------------------------------------------------------------------

function generate_io_equation_problem(ode::ODE, outputs)
    dim_x = length(ode.x_vars)
    if !isa(outputs, Array)
        outputs = [outputs]
    end

    # Creating a ring
    old_vars = map(string, gens(ode.poly_ring))
    var_names = vcat(
        old_vars,
        ["$(x)_dot" for x in ode.x_vars],
        ["y$(j)_$i" for i in 0:dim_x for j in 1:length(outputs)],
        ["$(u)_$i" for i in 1:dim_x for u in ode.u_vars],
    )
    ring, ring_vars = PolynomialRing(base_ring(ode.poly_ring), var_names)

    # Definiting a (partial) derivation on it
    derivation = Dict()
    for x in ode.x_vars
        derivation[str_to_var(string(x), ring)] = str_to_var(string(x) * "_dot", ring)
    end
    for i in 0:(dim_x - 1)
        for j in 1:length(outputs)
            derivation[str_to_var("y$(j)_$i", ring)] = str_to_var("y$(j)_$(i + 1)", ring)
        end
    end
    for u in ode.u_vars
        derivation[str_to_var(string(u), ring)] = str_to_var("$(u)_1", ring)
        for i in 1:(dim_x - 1)
            derivation[str_to_var("$(u)_$i", ring)] = str_to_var("$(u)_$(i + 1)", ring)
        end
    end

    # Generate equations
    x_equations = Dict()
    for x in ode.x_vars
        x_lifted = parent_ring_change(x, ring)
        num, den = map(p -> parent_ring_change(ode.poly_ring(p), ring), unpack_fraction(ode.equations[x]))
        x_equations[x_lifted] = den * derivation[x_lifted] - num
    end
    y_equations = []
    for (i, g) in enumerate(outputs)
        g_num, g_den = map(p -> parent_ring_change(ode.poly_ring(p), ring), unpack_fraction(g))
        push!(y_equations, g_den * str_to_var("y$(i)_0", ring) - g_num)
    end

    generic_point_generator = ODEPointGenerator(ode, outputs, ring)

    return (ring, derivation, x_equations, y_equations, generic_point_generator)
end

#------------------------------------------------------------------------------

function find_ioequation(ode::ODE, output, auto_var_change = true)
    """
    Find the io_equation of an ODE system
    Input:
        - ode, the ODE system
        - auto_var_change::Bool, whether or not to perform automatic variable change
    """
    #Initialization
    ring, derivation, x_equations, y_equations, point_generator = generate_io_equation_problem(ode, output)
    x_left = Set(keys(x_equations))
    y_orders = [0 for _ in y_equations]
 
    while true        
        var_degs = [(i, [degree(eq, x) for x in x_left if degree(eq, x) > 0]) for (i, eq) in enumerate(y_equations)]
        filter!(d -> length(d[2]) > 0, var_degs)
        if isempty(var_degs)
            return y_equations
        end
        @debug "Current degrees of io-equations $var_degs"
        @debug "Orders: $y_orders"
        @debug "Sizes: $([length(eq) for eq in y_equations])"

        # choosing the output to prolong
        outputs_with_scores = [
            (
                min(d[2]...),
                -count(x -> x == min(d[2]...), d[2]),
                length(y_equations[d[1]]),
                d[1]
            ) for d in var_degs               
        ]
        @debug "Scores: $outputs_with_scores"
        y_ind = sort(outputs_with_scores)[1][end]
        y_orders[y_ind] += 1
        @debug "Prolonging output number $y_ind"
        flush(stdout)

        #Calculate the Lie derivative of the io_relation
        @debug "Prolonging"
        next_y_equation = diff_poly(y_equations[y_ind], derivation)
        for x in x_left
            @debug "Eliminating the derivative of $x"
            next_y_equation = eliminate_var(x_equations[x], next_y_equation, derivation[x], point_generator)
        end
        
        #Choose variable to eliminate
        var_degs_next = [(degree(y_equations[y_ind], x), degree(next_y_equation, x), x) for x in x_left if degree(y_equations[y_ind], x) > 0]
        our_choice = sort(var_degs_next)[1]
        var_elim_deg, var_elim = our_choice[1], our_choice[3]
        
        @debug "Elimination of $var_elim, $(length(x_left) - 1) left; $(Dates.now())"
        flush(stdout)
        
        #Possible variable change for Axy + Bx + p(y) (x = var_elim)
        if auto_var_change && (var_elim_deg == 1)
            Ay_plus_B = coeff(y_equations[y_ind], [var_elim], [1])
            for x in setdiff(x_left, [var_elim])
                if degree(Ay_plus_B, x) == 1                      
                    A, B = divrem(Ay_plus_B, x)
                    A, B = simplify_frac(A, B)
                    if isempty(filter!(v -> (v in x_left), vars(A))) && (B != 0) #&& (length(coeffs(A))==1) 
                        #variable change x_i' --> x_i' - (B/A)', x_i --> x_i - B/A
                        @debug "\t Applying variable change: $(x) --> $(x) - ( $B )/( $A ); $(Dates.now())"
                        flush(stdout)
                        dB = diff_poly(B, derivation)
                        dA = diff_poly(A, derivation)
                        numer_d, denom_d = simplify_frac(A * dB - dA * B, A * A)
                        #!!!push x_i_dot first, otherwise there will be problem with plugin!!!
                        point_generator = generator_var_change(
                            point_generator, 
                            derivation[x], 
                            derivation[x] * denom_d + numer_d, denom_d
                        )
                        point_generator = generator_var_change(point_generator, x, A * x + B, A)

                        #Change current system
                        @debug "Change in the system"
                        x_equations[x] = make_substitution(x_equations[x], derivation[x], denom_d * derivation[x] - numer_d, denom_d)
                        for xx in x_left
                            @debug "\t Change in the equation for $xx"
                            x_equations[xx] = make_substitution(x_equations[xx], x, A * x - B, A)
                        end
                        @debug "Change in the outputs"
                        for i in 1:length(y_equations)
                            @debug "\t Change in the $i th output"
                            y_equations[i] = make_substitution(y_equations[i], x, A * x - B, A)
                        end
                        @debug "\t Change in the prolonged equation"
                        next_y_equation = make_substitution(next_y_equation, x, A * x - B, A)
                        #recalibrate system
                        @debug "Unmixing the derivatives"
                        for xx in setdiff(x_left, [x])
                            @debug "\t Unmixing $xx"
                            x_equations[x] = eliminate_var(x_equations[x], x_equations[xx], derivation[xx], point_generator)
                        end
                        @debug "Change of variables performed"
                        break
                    end
                end
            end  
        end

        #Eliminate var_elim from the system
        delete!(x_equations, var_elim)
        delete!(x_left, var_elim)
        @debug "Elimination in states"
        for x in x_left
            @debug "\t Elimination in the equation for $x"
            x_equations[x] = eliminate_var(x_equations[x], y_equations[y_ind], var_elim, point_generator)
        end
        @debug "Elimination in y_equations"
        for i in 1:length(y_equations)
            if i != y_ind
                @debug "Elimination in the $i th output"
                y_equations[i] = eliminate_var(y_equations[i], y_equations[y_ind], var_elim, point_generator)
            end
        end
        @debug "\t Elimination in the prolonged equation"
        y_equations[y_ind] = eliminate_var(y_equations[y_ind], next_y_equation, var_elim, point_generator)
    end
end
