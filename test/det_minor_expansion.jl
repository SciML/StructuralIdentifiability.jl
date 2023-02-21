@testset "Determinant by minor expansion" begin
    for d in 1:5
        for testcase in 1:10
            mat_space = Nemo.MatrixSpace(Nemo.QQ, d, d)
            rnd_matrix = mat_space([mod(rand(Int), 1000) for i in 1:d, j in 1:d])
            @test det(rnd_matrix) == det_minor_expansion(rnd_matrix)
        end
    end
end
