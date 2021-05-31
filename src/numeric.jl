"""
    Computes enough prolongations by Lie derivatives to use SIAN
    Current assumption: no inputs, one output, all polynomial
"""
function get_poly_system_Lie(ode::ODE{P}) where P <: MPolyElem{<: FieldElem}
    prolong_count = length(ode.x_vars) + length(ode.parameters)

    # creating differential ring
    R = DiffPolyRing(
        Nemo.QQ, 
        map(var_to_str, vcat(ode.x_vars, ode.y_vars)), 
        map(var_to_str, ode.parameters),
        vcat([0 for x in ode.x_vars], [prolong_count for y in ode.y_vars])
    )
    set_custom_derivations!(R, ode.x_equations)

    # equations forming square subsystem
    result_core = Array{P, 1}()
    # extra equations
    result_extra = Array{P, 1}()

    for y in ode.y_vars
        eq = y - ode.y_equations[y]
        push!(result_core, to_diffpoly(eq, R))
        for _ in 1:(prolong_count - 1)
            push!(result_core, diff(result_core[end], R))
        end
        push!(result_extra, diff(result_core[end], R))
    end

    return Dict(
        :diff_ring => R,
        :square_system => result_core, 
        :extra_equations => result_extra
    )
end
