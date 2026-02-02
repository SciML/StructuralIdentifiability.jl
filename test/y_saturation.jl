if GROUP == "All" || GROUP == "Core"
    @testset "Output saturation" begin
        # full nfkb model from benchmarks
        nfkb = @ODEmodel(
            x1'(t) = k_prod - k_deg * x1(t) - k1 * x1(t) * u(t),
            x2'(t) =
                -k3 * x2(t) - k_deg * x2(t) - a2 * x2(t) * x10(t) + t1 * x4(t) -
                a3 * x2(t) * x13(t) +
                t2 * x5(t) +
                (k1 * x1(t) - k2 * x2(t) * x8(t)) * u(t),
            x3'(t) = k3 * x2(t) - k_deg * x3(t) + k2 * x2(t) * x8(t) * u(t),
            x4'(t) = a2 * x2(t) * x10(t) - t1 * x4(t),
            x5'(t) = a3 * x2(t) * x13(t) - t2 * x5(t),
            x6'(t) = c6a * x13(t) - a1 * x6(t) * x10(t) + t2 * x5(t) - i1 * x6(t),
            x7'(t) = i1 * kv * x6(t) - a1 * x11(t) * x7(t),
            x8'(t) = c4 * x9(t) - c5 * x8(t),
            x9'(t) = c2 - c1 * x7(t) - c3 * x9(t),
            x10'(t) =
                -a2 * x2(t) * x10(t) - a1 * x10(t) * x6(t) + c4a * x12(t) - c5a * x10(t) -
                i1a * x10(t) + e1a * x11(t),
            x11'(t) = -a1 * x11(t) * x7(t) + i1a * kv * x10(t) - e1a * kv * x11(t),
            x12'(t) = c2a + c1a * x7(t) - c3a * x12(t),
            x13'(t) = a1 * x10(t) * x6(t) - c6a * x13(t) - a3 * x2(t) * x13(t) + e2a * x14(t),
            x14'(t) = a1 * x11(t) * x7(t) - e2a * kv * x14(t),
            x15'(t) = c2c + c1c * x7(t) - c3c * x15(t),
            y1(t) = x7(t),
            y2(t) = x10(t) + x13(t),
            y3(t) = x9(t),
            y4(t) = x1(t) + x2(t) + x3(t),
            y5(t) = x2(t),
            y6(t) = x12(t)
        )
        
        nfkb_orders = propose_orders(nfkb)
        new_nfkb = saturate_outputs(nfkb, nfkb_orders)
        @test length(new_nfkb.y_vars) == 14
        nfkb_time = @elapsed nfkb_ident_results = assess_identifiability(new_nfkb)
        @test nfkb_time < 500
        @test count(x -> x == :globally, values(nfkb_ident_results)) == 37
        @test count(x -> x == :nonidentifiable, values(nfkb_ident_results)) == 7

        # TumorPilis model from the benchmarks
        tum_pil = @ODEmodel(
                    T'(t) =
                        a * T(t) * (1 - b * T(t)) - c1 * N(t) * T(t) - D(t) * T(t) -
                        KT * M(t) * T(t), #tumor cells
                    L'(t) =
                        -m * L(t) - q * L(t) * T(t) - ucte * L(t)^2 +
                        r2 * C(t) * T(t) +
                        pI * L(t) * I(t) / (gI + I(t)) +
                        u1(t) - KL * M(t) * L(t), # tumor-specific effector cells, T-celss
                    N'(t) =
                        alpha1 - f * N(t) + g * T(t) * N(t) / (h + T(t)) - p * N(t) * T(t) - KN * M(t) * N(t), # non-specific effector cells, NK cells
                    C'(t) = alpha2 - beta * C(t) - KC * M(t) * C(t), #circulating lymphocytes
                    I'(t) =
                        pt * T(t) * L(t) / (gt + T(t)) + w * L(t) * I(t) - muI * I(t) + u2(t), # IL-2, VI = u2 aplicación directa, terapia de IL2
                    M'(t) = -gamma * M(t) + u1(t), #chemotherapy drug, terapia/aplicación de quimio, u1 = VM
                    y1(t) = L(t),
                    y2(t) = N(t),
                    y3(t) = M(t)
                )
        
        tum_pil_orders = propose_orders(tum_pil)
        new_tum_pil = saturate_outputs(tum_pil, tum_pil_orders)
        @test length(new_tum_pil.y_vars) == 6
        tum_pil_time = @elapsed tum_pil_results = assess_identifiability(new_tum_pil)
        @test tum_pil_time < 500
        @test count(x -> x == :globally, values(tum_pil_results)) == 22
        @test count(x -> x == :nonidentifiable, values(tum_pil_results)) == 9

    end
end
