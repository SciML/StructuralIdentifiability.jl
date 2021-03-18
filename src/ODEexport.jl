"""
    print_for_SIAN(ode)

Prints the ODE in the format accepted by SIAN (https://github.com/pogudingleb/SIAN)
"""
function print_for_SIAN(ode::ODE)
    varstr = Dict(x => var_to_str(x) * "(t)" for x in vcat(ode.x_vars, ode.u_vars, ode.y_vars))
    merge!(varstr, Dict(p => var_to_str(p) for p in ode.parameters))
    R_print, vars_print = PolynomialRing(base_ring(ode.poly_ring), [varstr[v] for v in gens(ode.poly_ring)])
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


