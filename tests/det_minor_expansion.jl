@testset "Determinant by minor expansion" begin
    for d = 1:5
        for testcase = 1:10
            mat_space = MatrixSpace(Nemo.QQ, d, d)
            rnd_matrix = mat_space([mod(rand(Int), 1000) for i = 1:d, j = 1:d])
            @test det(rnd_matrix) == det_minor_expansion(rnd_matrix)
        end
    end
end
