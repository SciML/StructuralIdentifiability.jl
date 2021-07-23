@testset "Power series matrix inverse" begin

    T, t = Nemo.PowerSeriesRing(Nemo.GF(2^31 - 1), 50, "t"; model=:capped_absolute)
    
    for d in 1:5
        S = Nemo.MatrixSpace(T, d, d)
        for case in 1:20
            M = S([random_ps(T) for i in 1:d, j in 1:d])
            while isequal(StructuralIdentifiability.LinearAlgebra.det(M), 0)
                M = S([random_ps(T) for i in 1:d, j in 1:d])
            end 
            invM = ps_matrix_inv(M)
            prod = invM * M
            @test prod == one(S)
        end
    end
end
