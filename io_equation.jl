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

function generate_io_equation_problem(ode::ODE, output)
    dim_x = length(ode.x_vars)

    # Creating a ring
    old_vars = map(string, gens(ode.poly_ring))
    var_names = vcat(
        old_vars,
        ["$(x)_dot" for x in ode.x_vars],
        ["y_$i" for i in 0:dim_x],
        ["$(u)_$i" for i in 1:dim_x for u in ode.u_vars],
    )
    ring, ring_vars = PolynomialRing(base_ring(ode.poly_ring), var_names)

    # Definiting a (partial) derivation on it
    derivation = Dict()
    for x in ode.x_vars
        derivation[str_to_var(string(x), ring)] = str_to_var(string(x) * "_dot", ring)
    end
    for i in 0:(dim_x - 1)
        derivation[str_to_var("y_$i", ring)] = str_to_var("y_$(i + 1)", ring)
    end
    for u in ode.u_vars
        derivation[str_to_var(string(u), ring)] = str_to_var("$(u)_1", ring)
        for i in 1:(dim_x - 1)
            derivation[str_to_var("$(u)_i", ring)] = str_to_var("$(u)_$(i + 1)", ring)
        end
    end

    # Generate equations
    x_equations = Dict()
    for x in ode.x_vars
        x_lifted = parent_ring_change(x, ring)
        num, den = map(p -> parent_ring_change(ode.poly_ring(p), ring), unpack_fraction(ode.equations[x]))
        x_equations[x_lifted] = den * derivation[x_lifted] - num
    end
    y_num, y_den = map(p -> parent_ring_change(ode.poly_ring(p), ring), unpack_fraction(output))
    y_equation = y_den * str_to_var("y_0", ring) - y_num

    # Construct generic point generator
    Lie_derivation = copy(derivation)
    for (x, f) in ode.equations
        f_num, f_den = map(p -> parent_ring_change(ode.poly_ring(p), ring), unpack_fraction(f))
        Lie_derivation[str_to_var("$x", ring)] = f_num // f_den
    end
    @debug "\t Computing Lie derivatives $(Dates.now())"
    flush(stdout)
    Lie_derivatives = []
    push!(Lie_derivatives, y_equation)
    for i in 1:dim_x
        push!(
            Lie_derivatives,
            unpack_fraction(diff_poly(Lie_derivatives[end], Lie_derivation))[1]
        )
    end
    generic_point_generator = RationalVarietyPointGenerator(
        vcat(collect(values(x_equations)), Lie_derivatives),
        map(s -> str_to_var(s, ring), vcat(old_vars, ["$(u)_$i" for i in 1:dim_x for u in ode.u_vars]))
    )

    return (ring, derivation, x_equations, y_equation, generic_point_generator)
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
    ring, derivation, x_equations, y_equation, point_generator = generate_io_equation_problem(ode, output)
    x_left = Set(keys(x_equations))
    order = 0
 
    while true        
        var_degs = [(x, degree(y_equation, x)) for x in x_left]
        filter!(d -> (d[2] > 0), var_degs)
        if isempty(var_degs)
            return y_equation
        end
        
        #Calculate the Lie derivative of the io_relation
        next_y_equation = diff_poly(y_equation, derivation)
        order += 1
        for x in x_left
            next_y_equation = eliminate_var(x_equations[x], next_y_equation, derivation[x], point_generator)
        end
        
        #Choose variable to eliminate
        var_degs_next = [(d[2], degree(next_y_equation, d[1]), d[1]) for d in var_degs]
        our_choice = sort(var_degs_next)[1]
        var_elim_deg, var_elim = our_choice[1], our_choice[3]
        
        @debug "Elimination of $var_elim, $(length(x_left)) left; $(Dates.now())"
        flush(stdout)
        
        #Possible variable change for Axy + Bx + p(y) (x = var_elim)
        if auto_var_change && (var_elim_deg == 1)
            Ay_plus_B = coeff(y_equation, [var_elim], [1])
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
                        x_equations[x] = make_substitution(x_equations[x], derivation[x], denom_d * derivation[x] - numer_d, denom_d)
                        for xx in x_left
                            x_equations[xx] = make_substitution(x_equations[xx], x, A * x - B, A)
                        end
                        y_equation = make_substitution(y_equation, x, A * x - B, A)
                        next_y_equation = make_substitution(next_y_equation, x, A * x - B, A)
                        #recalibrate system
                        for xx in setdiff(x_left, [x])
                            x_equations[x] = eliminate_var(x_equations[x], x_equations[xx], derivation[xx], point_generator)
                        end                        
                        break
                    end
                end
            end  
        end

        #Eliminate var_elim from the system
        delete!(x_equations, var_elim)
        delete!(x_left, var_elim)
        for x in x_left
            x_equations[x] = eliminate_var(x_equations[x], y_equation, var_elim, point_generator)
        end
        #Update io_relation
        y_equation = eliminate_var(y_equation, next_y_equation, var_elim, point_generator)
    end
    return y_equation
end
