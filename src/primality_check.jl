#------------------------------------------------------------------------------

function check_primality_zerodim(J::Singular.sideal{Singular.spoly{Singular.n_Q}})
    J = Singular.std(J)
    basis = gens(Singular.kbase(J))
    dim = length(basis)
    S = Nemo.MatrixSpace(Nemo.QQ, dim, dim)
    matrices = []

    for v in gens(base_ring(J))
        M = zero(S)
        for (i, vec) in enumerate(basis)
            image = Singular.reduce(v * vec, J)
            for (j, base_vec) in enumerate(basis)
                M[i, j] = Nemo.QQ(coeff(image, base_vec))
            end
        end
        push!(matrices, M)
    end
    generic_multiplication = sum([rand(1:100) * M for M in matrices])
    return isirreducible(Nemo.charpoly(generic_multiplication))
end

#------------------------------------------------------------------------------

function check_primality(polys::Dict{fmpq_mpoly, fmpq_mpoly}, extra_relations::Array{fmpq_mpoly, 1})
    leaders = collect(keys(polys))
    ring = parent(leaders[1])

    Rsing, vsing = Singular.PolynomialRing(Singular.QQ, [var_to_str(l) for l in leaders])
    eval_point = [v in keys(polys) ? v : ring(rand(1:100)) for v in gens(ring)]
    all_polys = vcat(collect(values(polys)), extra_relations)
    zerodim_ideal = Singular.Ideal(
        Rsing, 
        map(p -> parent_ring_change(AbstractAlgebra.evaluate(p, eval_point), Rsing), all_polys)
    )
    
    return check_primality_zerodim(zerodim_ideal)
end

#------------------------------------------------------------------------------

function check_primality(polys::Dict{fmpq_mpoly, fmpq_mpoly})
    return check_primality(polys, Array{fmpq_mpoly, 1}())
end

#------------------------------------------------------------------------------
