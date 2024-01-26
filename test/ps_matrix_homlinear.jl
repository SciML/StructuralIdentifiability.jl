@testset "Homogeneous linear differential equations" begin
    T, t =
        Nemo.power_series_ring(Nemo.Native.GF(2^31 - 1), 300, "t"; model = :capped_absolute)

    for d in 1:5
        for c in 1:5
            S = Nemo.matrix_space(T, d, d)
            A = random_ps_matrix(T, S)
            Sconst = Nemo.matrix_space(Nemo.Native.GF(2^31 - 1), d, d)
            Y0 = Sconst([rand(Int) % 100 for i in 1:d, j in 1:d])
            while StructuralIdentifiability.LinearAlgebra.det(Y0) == 0
                Y0 = Sconst([rand(Int) % 100 for i in 1:d, j in 1:d])
            end
            @time sol, sol_inv = ps_matrix_homlinear_de(A, Y0)
            to_be_zero = map(ps_diff, sol) - A * sol
            @test truncate_matrix(to_be_zero, 299) == truncate_matrix(zero(S), 299)
        end
    end
end
