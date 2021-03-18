@testset "Logarith of power series matrices" begin
    
    T, t = Nemo.PowerSeriesRing(Nemo.QQ, 300, "t"; model=:capped_absolute)
    
    for d in 1:5
        diag_elements = [1 + t * random_ps(T) for i in 1:d]
        S = Nemo.MatrixSpace(T, d, d)
        M = S([(i == j ? diag_elements[i] : T(0)) for i = 1:d, j = 1:d])
        result = ps_matrix_log(M)
        correct = S([(i == j ? log(diag_elements[i]) : T(0)) for i = 1:d, j = 1:d])
        @test result == correct
    end
    
    for c in 1:5
        f = random_ps(T) * t
        S = Nemo.MatrixSpace(T, 2, 2)
        m = S([
            T(1) f; 
            T(0) T(1)
        ])
        correct = S([
            T(0) f; 
            T(0) T(0)
        ])
        @test ps_matrix_log(m) == correct
    end
    
    for c in 1:5
        f, g = random_ps(T) * t, random_ps(T) * t
        S = Nemo.MatrixSpace(T, 3, 3)
        m = S([
            T(1) f    g; 
            T(0) T(1) f; 
            T(0) T(0) T(1)
        ])
        correct = S([
            T(0) f    (-1 // 2) * f^2 + g;
            T(0) T(0) f;
            T(0) T(0) T(0)
        ])
        @test ps_matrix_log(m) == correct
    end
    
end
