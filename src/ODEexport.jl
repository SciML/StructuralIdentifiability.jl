"""
    print_for_maple(ode, package)

Prints the ODE in the format accepted by maple packages
  - SIAN (https://github.com/pogudingleb/SIAN) if package=:SIAN
  - DifferentialAlgebra if package=:DifferentialAlgebra
  - DifferentialThomas if package=:DifferentialThomas
"""
function print_for_maple(ode::ODE, package=:SIAN)
    varstr = Dict(x => var_to_str(x) * "(t)" for x in vcat(ode.x_vars, ode.u_vars, ode.y_vars))
    merge!(varstr, Dict(p => var_to_str(p) for p in ode.parameters))
    R_print, vars_print = Nemo.PolynomialRing(base_ring(ode.poly_ring), [varstr[v] for v in gens(ode.poly_ring)])
    x_names = join(map(var_to_str, ode.x_vars), ", ")
    y_names = join(map(var_to_str, ode.y_vars), ", ")
    u_names = join(map(var_to_str, ode.u_vars), ", ")
    u_string = length(u_names) > 0 ? "[$u_names]" : ""
    ranking = "[$x_names], [$y_names], $u_string"

    result = ""

    # Forming the header
    if package == :SIAN
        result *= "read \"../IdentifiabilityODE.mpl\";\n\nsys := [\n"
    elseif package == :DifferentialAlgebra
        result *= "with(DifferentialAlgebra):\n"
        u_string = length(u_names) > 0 ? "[$u_names]" : ""
        result *= "ring_diff := DifferentialRing(blocks = [$ranking], derivations = [t]):\n"
        result *= "sys := [\n"
    elseif package == :DifferentialThomas
        result *= "with(DifferentialThomas):\nwith(Tools):\n"
        result *= "Ranking($ranking):\n"
        result *= "sys := [\n"
    else
        throw(Base.ArgumentError("Unknown package: $package"))
    end

    # Form the equations
    function _rhs_to_str(rhs)
        num, den = unpack_fraction(rhs)
        res = string(evaluate(num, vars_print))
        if den != 1
             res = "($res) / ($(evaluate(den, vars_print)))"
        end
        return res
    end

    eqs = []
    for (x, f) in ode.x_equations
        push!(eqs, ("diff(" * var_to_str(x) * "(t), t)", _rhs_to_str(f)))
    end
    for (y, g) in ode.y_equations
        push!(eqs, (var_to_str(y) * "(t)", _rhs_to_str(g)))
    end
    if package == :SIAN 
        result *= join(map(a -> a[1] * " = " * a[2], eqs), ",\n") * "\n];\n"
    else
        result *= join(map(a -> "$(a[1]) - ($(a[2]))", eqs), ",\n") * "\n];\n"
    end

    # Forming the footer
    if package == :SIAN
        result *= "CodeTools[CPUTime](IdentifiabilityODE(sys, GetParameters(sys)));"
    elseif package == :DifferentialAlgebra
        result *= "CodeTools[CPUTime](RosenfeldGroebner(sys, ring_diff));"
    else
        result *= "CodeTools[CPUTime](ThomasDecomposition(sys));"
    end

    # Eliminating Julia integer divisions and variable name I
    result = replace(result, "//" => "/")
    result = replace(result, "I(t)" => "II(t)")
    return result
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

    result = result * "FLAG_:=1\$\nSHOWTIME\$\nDAISY()\$\nSHOWTIME\$\nEND\$\n"

    result = replace(result, "//" => "/")
    return result
end
