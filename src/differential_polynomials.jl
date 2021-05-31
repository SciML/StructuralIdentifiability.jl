function diffvar(name::String, ord::Int)
    return "$(name)_der_$ord"
end

#------------------------------------------------------------------------------

struct DiffPolyRing
    diff_var_names::Array{String, 1},
    max_orders::Array{Int, 1},
    parameters::Array{<: MPolyElem, 1}
    derivation::Dict{<: MPolyElem, <: MPolyElem}
    ring::MPolyRing

    function DiffPolyRing(
        coef_field::Ring,
        diff_var_names::Array{String, 1}, 
        param_names::Array{String, 1}, 
        max_orders::Array{Int, 1}
    )
        varnames = param_names
        for (i, name) in enumerate(diff_var_names)
            for ord in 0:max_orders[i]
                push!(varnames, diffvar(name, ord))
            end
        end
        ring, _ = PolynomialRing(coef_field, varnames)
        parameters = map(p -> str_to_var(p, ring), param_names)
        derivation = Dict{MPolyElem, MPolyElem}()
        for (i, v) in enumerate(diff_var_names)
            for ord in 0:(max_orders[i] - 1)
                derivation[str_to_var(diffvar(v, ord), ring)] = str_to_var(diffvar(v, ord + 1), ring)
            end
        end

        return new(diff_var_names, max_orders, parameters, derivation, ring)
    end
end

#------------------------------------------------------------------------------

"""
    Converts a polynomial to a differential one by assigning zero derivatives to all the variables
"""
function to_diffpoly(p::MPolyElem, R::DiffPolyRing)
    eval_from = Array{MPolyElem, 1}()
    eval_to = Array{MPolyElem, 1}()
    for v in vars(p)
        push!(eval_from, v)

        s = var_to_str(v)
        if s in R.diff_var_names
            push!(eval_to, str_to_var(diffvar(s, 0), R.ring))
        elseif s in map(var_to_str, R.parameters)
            push!(eval_to, str_to_var(s, R.ring))
        else
            push!(eval_to, zero(R.ring))
        end
    end
    return evaluate(p, eval_from, eval_to)
end

#------------------------------------------------------------------------------

function diff(p::MPolyElem, R::DiffPolyRing)
    @assert parent(p) == R.ring

    result = zero(R.ring)
    for (v, d) in R.derivation
        result += derivative(p, v) * d
    end
    return result
end

#------------------------------------------------------------------------------

function set_custom_derivations!(R::DiffPolyRing, ders::Dict)
    for (v, d) in ders
        @assert var_to_str(v) in R.diff_var_names
        R.derivation[to_diffpoly(v, R)] = to_diffpoly(d, R)
    end
end
