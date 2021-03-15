@testset "Homogeneous linear differential equations" begin
    
    T, t = Nemo.PowerSeriesRing(Nemo.GF(2^31 - 1), 300, "t"; model=:capped_absolute)
    
    for d = 1:5
        for c = 1:5
            S = Nemo.MatrixSpace(T, d, d)
            A = random_ps_matrix(T, S)
            Sconst = Nemo.MatrixSpace(Nemo.GF(2^31 - 1), d, d)
            Y0 = Sconst([rand(Int) % 100 for i = 1:d, j = 1:d])
            @time sol, sol_inv = ps_matrix_homlinear_de(A, Y0)
            to_be_zero = map(ps_diff, sol) - A * sol
            @test truncate_matrix(to_be_zero, 299) == truncate_matrix(zero(S), 299)
        end
    end
    
end
