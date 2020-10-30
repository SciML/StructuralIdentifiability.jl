using Nemo

function check_primality_univariate(polys::Array{fmpq_poly, 1})
    degrees = map(degree, polys)
    basis = vec(collect(Base.Iterators.product([0:(d - 1) for d in degrees]...)))
    S = MatrixSpace(Nemo.QQ, length(basis), length(basis))
    matrices = []
    for i in 1:length(polys)
        d = degrees[i]
        M = zero(S)
        e = tuple([j == i ? 1 : 0 for j in 1:length(polys)]...)
        for (j, b) in enumerate(basis)
            if b[i] < d - 1
                M[findfirst(v -> v == b .+ e, basis), j] = 1
            else
                for k in 1:d
                    M[findfirst(v -> v == b, basis), j] = -coeff(polys[i], d - k) // coeff(polys[i], d)
                    b = b .- e
                end
            end
        end
        push!(matrices, M)
    end
    generic_multiplication = sum([rand(1:100) * M for M in matrices])
    return isirreducible(Nemo.charpoly(generic_multiplication))
end

function check_primality(polys::Dict{fmpq_mpoly, fmpq_mpoly})
    leaders = collect(keys(polys))
    ring = parent(leaders[1])
    R, x_aux = PolynomialRing(Nemo.QQ, "x")
    eval_point = [v in keys(polys) ? x_aux : R(rand(1:100)) for v in gens(ring)]
    univar = [evaluate(p, eval_point) for p in values(polys)]
    return check_primality_univariate(univar)
end
