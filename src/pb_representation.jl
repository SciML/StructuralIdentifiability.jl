"""
The structure for storing a projection-based representation of differential ideal
(see Section 2.3 https://arxiv.org/abs/2111.00991).
Contains the following fields:
-  `y_names` - the names of the variables with finite order in the profile (typically, outputs)
-  `u_names` - the names of the variables with infinite order in the profile (typically, inputs)
-  `param_names` - the names of the parameters
-  `profile` - the profile of the PB-representation (see Definiton 2.13) as a dict from `y_names` with finite orders to the orders
-  `projections` - the corresponding projections (see Definition 2.15) as a dict from `y_names` to the projections
"""
struct PBRepresentation
    y_names::Array{String} # variables with finite orders in the profile
    u_names::Array{String} # variables with infinite orders in the profile
    param_names::Array{String} # scalar parameters
    profile::Dict{String, Int} # the profile restricted on the y-variables
    projections::Dict{String, <:MPolyElem}

    function PBRepresentation(ode::ODE, io_equations)
        if length(keys(io_equations)) > length(ode.y_vars)
            throw(DomainError("This projection-based representation is not a charset, not yet implemented"))
        end
        y_names = map(var_to_str, ode.y_vars)
        u_names = map(var_to_str, ode.u_vars)
        param_names = map(var_to_str, ode.parameters)
        old_ring = parent(first(values(io_equations)))
        new_varnames = filter(
            v -> (v in param_names) || decompose_derivative(v, vcat(y_names, u_names)) != nothing,
            map(var_to_str, gens(old_ring))
        )
        newring, _ = Nemo.PolynomialRing(base_ring(old_ring), new_varnames)

        profile = Dict{String, Int}()
        projections = Dict{String, MPolyElem}()
        for (y, eq) in io_equations
            (name, ord) = decompose_derivative(var_to_str(y), y_names)
            profile[name] = ord
            projections[name] = parent_ring_change(eq, newring)
        end
        new(y_names, u_names, param_names, profile, projections)
    end

    # for testing
    function PBRepresentation(y_names, u_names, param_names, profile, projections)
        new(y_names, u_names, param_names, profile, projections)
    end
end

# -----------------------------------------------------------------------------

"""
    find_leader(vars, pbr)

Among the variables `vars`, determines the leading derivative if the y-variable
(if exists) with respect to the ordering defined by the PB-representation
(see Remark 2.20 in https://arxiv.org/abs/2111.00991)
"""
function find_leader(vars::Array{<: MPolyElem}, pbr::PBRepresentation)
    y_ders = filter(v -> decompose_derivative(var_to_str(v), pbr.y_names) != nothing, vars)
    if length(y_ders) == 0
        return nothing
    end
    y_ders_ext = [(decompose_derivative(var_to_str(y), pbr.y_names), y) for y in y_ders]
    return sort(
        y_ders_ext, rev = true, 
        by = p -> (
            p[1][2] - pbr.profile[p[1][1]], 
            findfirst(n -> n == p[1][1], pbr.y_names)
        )
    )[1][2]
end

# -----------------------------------------------------------------------------

"""
    common_ring(poly, pbr)

For a polynomail `poly` in the same differential variables as `pbr`, finds
a polynomial ring sufficient for carrying out the reduction and the 
corresponding differentiation mapping on the variables
"""
function common_ring(poly::MPolyElem, pbr::PBRepresentation)
    max_ords = Dict{String, Int}(v => 0 for v in vcat(pbr.y_names, pbr.u_names))
    new_params = Array{String, 1}()
    for v in vars(poly)
        d = decompose_derivative(var_to_str(v), vcat(pbr.y_names, pbr.u_names))
        if d != nothing
            max_ords[d[1]] = max(d[2], max_ords[d[1]])
        elseif !(var_to_str(v) in pbr.param_names)
            @warn "New variable $(var_to_str(v)), treating as a scalar parameter"
            push!(new_params, var_to_str(v))
        end
    end

    max_offset = max([max_ords[y] - pbr.profile[y] for y in pbr.y_names]...)

    varnames = Array{String, 1}()
    for y in pbr.y_names
        append!(varnames, ["$(y)_$h" for h in 0:(max_offset + pbr.profile[y])])
    end
    for u in pbr.u_names
        append!(
                varnames, 
                ["$(u)_$h" for h in 0:max(max_ords[u], 
                max_offset + max([difforder(p, u) for p in values(pbr.projections)]...))]
        )
    end
    append!(varnames, pbr.param_names)
    append!(varnames, new_params)

    newring, _ = StructuralIdentifiability.Nemo.PolynomialRing(base_ring(parent(poly)), varnames)
    derivation = Dict{MPolyElem, MPolyElem}()
    for v in varnames
        d = decompose_derivative(v, vcat(pbr.y_names, pbr.u_names))
        if d == nothing
            derivation[str_to_var(v, newring)] = zero(newring)
        else
            der = "$(d[1])_$(d[2] + 1)"
            if der in varnames
                derivation[str_to_var(v, newring)] = str_to_var(der, newring)
            else
                derivation[str_to_var(v, newring)] = zero(newring)
            end
        end
    end
    return (newring, derivation)
end

# -----------------------------------------------------------------------------

"""
    lc_univariate(f, x)

Computes the leading coefficient of `f` viewed as a univariate polynomail in variable `x`
"""
function lc_univariate(f::MPolyElem, x::MPolyElem)
    FieldType = typeof(one(base_ring(parent(f))))
    dict_result = Dict{Array{Int,1}, FieldType}()
    x_ind = findfirst(v -> v == x, gens(parent(f)))
    cur_deg = 0
    for (monom, coef) in zip(exponent_vectors(f), coefficients(f)) 
        if monom[x_ind] > cur_deg
            cur_deg = monom[x_ind]
            dict_result = Dict{Array{Int, 1}, FieldType}()
        end
        if monom[x_ind] == cur_deg
            monom[x_ind] = 0
            dict_result[monom] = coef
        end
    end
    return dict_to_poly(dict_result, parent(f))
end

# -----------------------------------------------------------------------------

"""
    pseudodivision(f, g, x)

Computes the result of pseudodivision of `f` by `g` as univariate polynomials in `x`
Input:
-  `f` - the polynomail to be divided
-  `g` - the polynomial to divide by
-  `x` - the variable for the division

Output: the pseudoreminder of `f` divided by `g` w.r.t. `x`
"""
function pseudodivision(f::MPolyElem, g::MPolyElem, x::MPolyElem)
    result = f
    lcg = lc_univariate(g, x)
    while Nemo.degree(result, x) >= Nemo.degree(g, x)
        degdiff = Nemo.degree(result, x) - Nemo.degree(g, x)
        result = lcg * result - lc_univariate(result, x) * x^degdiff * g
        lcgcd = gcd(result, lcg)
        if total_degree(lcgcd) != 0
            result = divexact(result, lcgcd)
        end
    end
    return result
end

# -----------------------------------------------------------------------------

function diff(p::MPolyElem, derivation::Dict{<: MPolyElem, <: MPolyElem}, i::Int)
    if i == 0
        return p
    end
    if i == 1
        return sum([derivative(p, x) * derivation[x] for x in vars(p)])
    end
    return diff(diff(p, derivation, 1), derivation, i - 1)
end

# -----------------------------------------------------------------------------

"""
    diffreduce(diffpoly, pbr)

Computes the result of differential reduction of a differential polynomial
`diffpoly` with respect to the charset defined by a PB-representation `pbr`
Input:
-  `diffpoly` - a polynomial representing a differential polynomial to be reduced
-  `pbr` - a projection-based representation

Output: the result of differential reduction of `diffpoly` by `pbr` considered as a characteristic set (see Remark 2.20 in the paper)
"""
function diffreduce(diffpoly::MPolyElem, pbr::PBRepresentation)
    (ring, der) = common_ring(diffpoly, pbr)

    result = parent_ring_change(diffpoly, ring)
    ext_projections = Dict(y => parent_ring_change(f, ring) for (y, f) in pbr.projections)

    while true
        lead = find_leader(vars(result), pbr)
        if lead == nothing
            return result
        end
        (var, ord) = decompose_derivative(var_to_str(lead), pbr.y_names)
        if ord < pbr.profile[var]
            return result
        end
        if ord == pbr.profile[var] && Nemo.degree(result, lead) < Nemo.degree(ext_projections[var], lead)
            return result
        end
        reducer = diff(ext_projections[var], der, ord - pbr.profile[var])
        result = pseudodivision(result, reducer, lead)
    end
end

# -----------------------------------------------------------------------------

"""
    io_switch(pbr)

In a single-output pb-representation `pbr` makes the leading variable to be the first of the inputs
"""
function io_switch!(pbr::PBRepresentation)
    if length(pbr.y_names) > 1
        throw(ArgumentError("Not implemented for multiple outputs"))
    end
    diffpoly = first(values(pbr.projections))
    u = first(pbr.u_names)
    ordu = -1
    for v in vars(diffpoly)
        d = decompose_derivative(var_to_str(v), [u])
        if d != nothing
            ordu = max(d[2], ordu)
        end
    end
    y = first(pbr.y_names)
    pbr.u_names[1] = y
    pbr.y_names[1] = u
    pbr.projections[u] = diffpoly
    delete!(pbr.projections, y)
    pbr.profile[u] = ordu
    delete!(pbr.profile, y)
end
