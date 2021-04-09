"""
    print_for_SIAN(ode)

Prints the ODE in the format accepted by SIAN (https://github.com/pogudingleb/SIAN)
"""
function print_for_SIAN(ode::ODE)
    varstr = Dict(x => var_to_str(x) * "(t)" for x in vcat(ode.x_vars, ode.u_vars, ode.y_vars))
    merge!(varstr, Dict(p => var_to_str(p) for p in ode.parameters))
    R_print, vars_print = Nemo.PolynomialRing(base_ring(ode.poly_ring), [varstr[v] for v in gens(ode.poly_ring)])
    result = "read '../IdentifiabilityODE.mpl';\n\nsigma := [\n"

    function _rhs_to_str(rhs)
        num, den = unpack_fraction(rhs)
        result = string(evaluate(num, vars_print))
        if den != 1
             result = "($result) / ($(evaluate(den, vars_print)))"
        end
        return result
    end

    eqs = []
    for (x, f) in ode.x_equations
        push!(eqs, "diff(" * var_to_str(x) * "(t), t) = $(_rhs_to_str(f))")
    end
    for (y, g) in ode.y_equations
        push!(eqs, var_to_str(y) * "(t) = $(_rhs_to_str(g))")
    end
    return result * join(eqs, ",\n") * "\n];\nIdentifiabilityODE(sigma, GetParameters(sigma));"
end

#------------------------------------------------------------------------------

"""
    print_for_DAISY(ode)

Prints the ODE in the format accepted by DAISY (https://daisy.dei.unipd.it/)
"""
function print_for_DAISY(ode::ODE)
    result = ""

    # listing functions (ordering u, y, x matters!)
    result = result * "B_:={" * join(map(var_to_str, vcat(ode.u_vars, ode.y_vars, ode.x_vars)), ", ") * "}\$\n"

    result = result * "FOR EACH EL_ IN B_ DO DEPEND EL_,T\$\n\n"

    # listing parameters
    result = result * "B1_:={" * join(map(var_to_str, ode.parameters), ", ") * "}\$\n"

    result = result * " %NUMBER OF STATES\nNX_:=$(length(ode.x_vars))\$\n"
    result = result * " %NUMBER OF INPUTS\nNU_:=$(length(ode.u_vars))\$\n"
    result = result * " %NUMBER OF OUTPUTS\nNY_:=$(length(ode.y_vars))\$\n\n"

    eqs = []

    function _lhs_to_str(lhs)
        num, den = unpack_fraction(lhs)
        rslt = string(num)
        if den != 1
             rslt = "($rslt) / ($den)"
        end
        return rslt
    end

    for (x, f) in ode.x_equations
        push!(eqs, "df($(var_to_str(x)), t) = " * _lhs_to_str(f))
    end
    for (y, g) in ode.y_equations
        push!(eqs, "$(var_to_str(y)) = " * _lhs_to_str(g))
    end

    result = result * "C_:={" * join(eqs, ",\n") * "}\$\n"

    result = result * "SEED_:=25\$\nDAISY()\$\nEND\$\n"

    return result
end
