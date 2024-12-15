# ------------------------------------------------------------------------------

function check_primality_zerodim(J::Array{QQMPolyRingElem, 1})
    J = Groebner.groebner(J)
    basis = Groebner.kbase(J)
    dim = length(basis)
    S = Nemo.matrix_space(Nemo.QQ, dim, dim)
    matrices = []
    @debug "$J $basis"
    @debug "Dim is $dim"
    for v in gens(parent(first(J)))
        M = zero(S)
        for (i, vec) in enumerate(basis)
            image = Groebner.normalform(J, v * vec)
            for (j, base_vec) in enumerate(basis)
                M[i, j] = Nemo.QQ(coeff(image, base_vec))
            end
        end
        push!(matrices, M)
        @debug "Multiplication by $v: $M"
    end
    generic_multiplication = sum(Nemo.QQ(rand(1:100)) * M for M in matrices)
    @debug generic_multiplication

    R, t = Nemo.polynomial_ring(Nemo.QQ, "t")
    @debug "$(Nemo.charpoly(R, generic_multiplication))"

    return Nemo.is_irreducible(Nemo.charpoly(R, generic_multiplication))
end

#------------------------------------------------------------------------------
"""
    check_primality(polys::Dict{QQMPolyRingElem, QQMPolyRingElem}, extra_relations::Array{QQMPolyRingElem, 1})

The function checks if the ideal generated by the polynomials and saturated at
the leading coefficient with respect to the corresponding variables is prime
over rationals.

The `extra_relations` allows adding more polynomials to the generators (not affecting the saturation).
"""
function check_primality(
    polys::Dict{QQMPolyRingElem, QQMPolyRingElem},
    extra_relations::Array{QQMPolyRingElem, 1},
)
    leaders = collect(keys(polys))
    ring = parent(leaders[1])

    Rspec, vspec = Nemo.polynomial_ring(Nemo.QQ, [var_to_str(l) for l in leaders])
    eval_point = [v in keys(polys) ? v : ring(rand(1:100)) for v in gens(ring)]
    all_polys = vcat(collect(values(polys)), extra_relations)
    zerodim_ideal =
        collect(map(p -> parent_ring_change(evaluate(p, eval_point), Rspec), all_polys))

    return check_primality_zerodim(zerodim_ideal)
end

#------------------------------------------------------------------------------
"""
    check_primality(polys::Dict{QQMPolyRingElem, QQMPolyRingElem})

The function checks if the ideal generated by the polynomials and saturated at
the leading coefficient with respect to the corresponding variables is prime
over rationals.
"""
function check_primality(polys::Dict{QQMPolyRingElem, QQMPolyRingElem})
    return check_primality(polys, Array{QQMPolyRingElem, 1}())
end

# ------------------------------------------------------------------------------
