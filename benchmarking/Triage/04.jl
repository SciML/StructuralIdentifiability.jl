# https://doi.org/10.1016/j.automatica.2012.09.015
# alpha := 1

using StructuralIdentifiability, Nemo

# Eq. (9), 2m states
function make_spr9(m::Integer)
    names = vcat(
        ["C$i(t)" for i in 1:m],
        ["B$i(t)" for i in 1:m],
        ["ka$i" for i in 1:m],
        ["kd$i" for i in 1:m],
        ["CT$i" for i in 1:m],
        ["R", "kM", "h", "u(t)", "y(t)"]
    )

    S, vars = polynomial_ring(QQ, names)
    v = Dict(zip(names, vars))

    C = [v["C$i(t)"] for i in 1:m]
    B = [v["B$i(t)"] for i in 1:m]
    ka = [v["ka$i"] for i in 1:m]
    kd = [v["kd$i"] for i in 1:m]
    CT = [v["CT$i"] for i in 1:m]

    R, kM, h = v["R"], v["kM"], v["h"]
    u, y = v["u(t)"], v["y(t)"]

    A = R - sum(B)

    P = typeof(vars[1])
    F = typeof(vars[1] // one(S))
    xeqs = Dict{P, F}()

    for i in 1:m
        reaction = ka[i]*C[i]*A - kd[i]*B[i]

        xeqs[C[i]] = (
            -reaction + kM*(CT[i]*u - C[i])
        ) // h

        xeqs[B[i]] = reaction // one(S)
    end

    yeqs = Dict{P, F}(y => sum(B) // one(S))

    return StructuralIdentifiability.ODE{P}(vcat(C, B), [y], xeqs, yeqs, [u])
end

# Eq. (12), m states
function make_spr12(m::Integer)
    names = vcat(
        ["B$i(t)" for i in 1:m],
        ["ka$i" for i in 1:m],
        ["kd$i" for i in 1:m],
        ["CT$i" for i in 1:m],
        ["R", "kM", "u(t)", "y(t)"]
    )

    S, vars = polynomial_ring(QQ, names)
    v = Dict(zip(names, vars))

    B = [v["B$i(t)"] for i in 1:m]
    ka = [v["ka$i"] for i in 1:m]
    kd = [v["kd$i"] for i in 1:m]
    CT = [v["CT$i"] for i in 1:m]

    R, kM = v["R"], v["kM"]
    u, y = v["u(t)"], v["y(t)"]

    A = R - sum(B)

    P = typeof(vars[1])
    F = typeof(vars[1] // one(S))
    xeqs = Dict{P, F}()

    for i in 1:m
        xeqs[B[i]] = (
            kM*(ka[i]*CT[i]*u*A - kd[i]*B[i])
        ) // (kM + ka[i]*A)
    end

    yeqs = Dict{P, F}(y => sum(B) // one(S))

    return StructuralIdentifiability.ODE{P}(B, [y], xeqs, yeqs, [u])
end

spr9 = make_spr9(2)
assess_identifiability(spr9)

spr12 = make_spr12(2)
assess_identifiability(spr12)
