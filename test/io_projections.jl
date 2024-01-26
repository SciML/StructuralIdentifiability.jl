@testset "IO-projections (+ extra projection)" begin
    cases = []

    # Example from remark 3 in https://arxiv.org/pdf/2111.00991.pdf
    ode = @ODEmodel(
        x1'(t) = (1 + x1(t)^2) // 2,
        x2'(t) = (1 - x1(t)^2) // (1 + x1(t)^2),
        y1(t) = 2 * x1(t) // (b * (1 + x1(t)^2)),
        y2(t) = x2(t)
    )
    push!(cases, ode)

    #--------------------------------------------------------

    # Example from https://github.com/SciML/StructuralIdentifiability.jl/issues/132
    ode = @ODEmodel(
        Complex'(t) =
            1 / C * (
                ((2 * kon * free_receptor(t) * Drug(t) - koff * Complex(t)) * C) -
                (ke_Complex * Complex(t) * C) - (
                    (kon_2 * free_receptor(t) * Complex(t) - 2 * koff * Complex_2(t)) * C
                )
            ),
        Complex_2'(t) =
            1 / C * (
                (
                    (kon_2 * free_receptor(t) * Complex(t) - 2 * koff * Complex_2(t)) * C
                ) - (ke_Complex_2 * Complex_2(t) * C)
            ),
        Drug'(t) =
            1 / C * (
                -(ke_Drug * Drug(t) * C) -
                ((2 * kon * free_receptor(t) * Drug(t) - koff * Complex(t)) * C) +
                (45 / 100 * ka * Drug_SC(t))
            ),
        free_receptor'(t) =
            1 / C * (
                -((2 * kon * free_receptor(t) * Drug(t) - koff * Complex(t)) * C) +
                (66 / 2500 * C) - (kdeg_free_receptor * free_receptor(t) * C) - (
                    (kon_2 * free_receptor(t) * Complex(t) - 2 * koff * Complex_2(t)) * C
                )
            ),
        Drug_SC'(t) =
            -(45 / 100 * ka * Drug_SC(t)) - ((1 - 45 / 100) * ka * Drug_SC(t)) +
            u_SC(t),
        y1(t) = Drug(t),
        y2(t) = free_receptor(t) + Complex(t) + 2 * Complex_2(t)
    )
    push!(cases, ode)

    #---------------------------------------------------------

    for ode in cases
        proj, gpg, _ = find_ioprojections(ode, false)
        for p in values(proj)
            @test choose([p], gpg) == p
        end
        @test !check_primality(proj)
        # taking simply a sum instead of random linear combination
        extra_projection = sum(keys(proj))
        proj, gpg, projection_poly = find_ioprojections(ode, false, extra_projection)
        @test choose([projection_poly], gpg) == projection_poly
        @test check_primality(proj, [projection_poly])
    end
end
