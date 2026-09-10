include(joinpath(@__DIR__, "..", "shared", "test_setup.jl"))

if GROUP == "All" || GROUP == "Core"
    @testset "Primality check (zerodim subroutine)" begin
        # check_primality_zerodim goes through Nemo finite-field matrices that
        # fail to construct on i686 (InexactError / FqField) until Nemocas/Nemo#2358.
        if Sys.WORD_SIZE != 64
            @info "Skipping primality zerodim tests on $(Sys.WORD_SIZE)-bit (Nemo#2358)"
            return
        end

        R, (x, y) = Nemo.polynomial_ring(Nemo.QQ, ["x", "y"])

        @test check_primality_zerodim([x^2 - 1, y^2 - 4]) == false

        @test check_primality_zerodim([(x + 5) * (x^3 - 7), y - 3]) == false

        @test check_primality_zerodim([x^3 - 5, y - 1]) == true

        @test check_primality_zerodim([x^2 + 1, y^3 - 3 * x + x + 5]) == true

        # not prime over any modulous but prime over Q
        @test check_primality_zerodim([x, y^4 + 1]) == true
    end
end
