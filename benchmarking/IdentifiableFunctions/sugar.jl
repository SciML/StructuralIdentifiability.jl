using Nemo

function spolynomial(f, g)
    mf, mg = leading_monomial(f), leading_monomial(g)
    m = lcm(mf, mg)
    co_mf, co_mg = divexact(m, mf), divexact(m, mg)
    cf, cg = leading_coefficient(f), leading_coefficient(g)
    c = gcd(cf, cg)
    co_cf, co_cg = divexact(cg, c), divexact(cf, c)
    spoly = co_cf * co_mf * f - co_cg * co_mg * g
    spoly
end

function buchberger_criterion(basis, critical_pairs, i, j)
    fi, fj = basis[i], basis[j]
    m = gcd(leading_monomial(fi), leading_monomial(fj))
    isone(m)
end

ext(i, j) = (min(i, j), max(i, j))
function staircase_criterion(basis, critical_pairs, i, j)
    fi, fj = basis[i], basis[j]
    m = lcm(leading_monomial(fi), leading_monomial(fj))
    for k in eachindex(basis)
        k == i || k == j && continue
        !(ext(i, k) in critical_pairs) && !(ext(j, k) in critical_pairs) && continue
        if first(divides(m, leading_monomial(basis[k])))
            return true
        end
    end
    false
end

lcm_deg(i, j, basis) = total_degree(lcm(basis[i], basis[j]))
function sugar_deg(i, j, basis, sugar_cubes)
    f, g = basis[i], basis[j]
    mf, mg = leading_monomial(f), leading_monomial(g)
    m = lcm(mf, mg)
    co_mf, co_mg = divexact(m, mf), divexact(m, mg)
    max(total_degree(co_mf) + sugar_cubes[i], total_degree(co_mg) + sugar_cubes[j])
end

function select_critical_pair!(basis, critical_pairs, sugar_cubes, use_sugar)
    if use_sugar
        _, j = findmin(x -> sugar_deg(x[1], x[2], basis, sugar_cubes), critical_pairs)
    else
        _, j = findmin(x -> lcm_deg(x[1], x[2], basis), critical_pairs)
    end
    pair = critical_pairs[j]
    deleteat!(critical_pairs, j)
    deg =
        use_sugar ? sugar_deg(pair[1], pair[2], basis, sugar_cubes) :
        lcm_deg(pair[1], pair[2], basis)
    pair, deg
end

function update!(critical_pairs, basis, sugar_cubes, j)
    for i in 1:(j - 1)
        if !buchberger_criterion(basis, critical_pairs, i, j)
            push!(critical_pairs, (i, j))
        end
    end
    nothing
end

function normalform(p, basis, sugar_cubes)
    q, r = divrem(p, basis)
    @assert length(basis) == length(sugar_cubes)
    sugar = maximum(map(i -> total_degree(q[i]) + sugar_cubes[i], 1:length(sugar_cubes)))
    r, sugar
end

function gb(polys; use_sugar=true)
    # The classic Buchberger algorithm
    sort!(polys, by=leading_monomial)
    critical_pairs = Vector{NTuple{2, Int}}()
    basis = empty(polys)
    sugar_cubes = Int[]
    for poly in polys
        push!(basis, poly)
        push!(sugar_cubes, total_degree(poly))
        update!(critical_pairs, basis, sugar_cubes, length(basis))
    end
    d = 0
    while !isempty(critical_pairs)
        ((i, j), sugar1) =
            select_critical_pair!(basis, critical_pairs, sugar_cubes, use_sugar)
        buchberger_criterion(basis, critical_pairs, i, j) && continue
        # staircase_criterion(basis, critical_pairs, i, j) && continue
        if iszero(d % 10)
            deg = if use_sugar
                sugar_deg(i, j, basis, sugar_cubes)
            else
                lcm_deg(i, j, basis)
            end
            @info "Critical pair of degree $(deg), pairs left: $(length(critical_pairs))"
        end
        spoly = spolynomial(basis[i], basis[j])
        nf, sugar2 = normalform(spoly, basis, sugar_cubes)
        sugar = max(sugar1, sugar2)
        iszero(nf) && continue
        push!(basis, nf)
        push!(sugar_cubes, sugar)
        update!(critical_pairs, basis, sugar_cubes, length(basis))
        d += 1
    end
    basis
end
