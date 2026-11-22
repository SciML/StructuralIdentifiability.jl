include(joinpath(@__DIR__, "..", "shared", "test_setup.jl"))

if GROUP == "All" || GROUP == "Core"
    @testset "Determinant by minor expansion" begin
        for d in 1:5
            for testcase in 1:10
                mat_space = Nemo.matrix_space(Nemo.QQ, d, d)
                rnd_matrix = mat_space([mod(rand(Int), 1000) for i in 1:d, j in 1:d])
                @test det(rnd_matrix) == det_minor_expansion(rnd_matrix)
            end
        end
    end

    @testset "Common polynomial factor in a determinant row" begin
        R, vars = Nemo.polynomial_ring(Nemo.QQ, ["x$i" for i in 1:8])
        x, y, z = vars[1:3]
        shared = x^2 * y^2 * (x + y)^3
        a, b = 6 * y^3, x^3 * (2 * x + 3 * y)
        entries = Nemo.matrix(R, [a * shared b * shared])

        _, factors = StructuralIdentifiability.simplify_matrix(entries)
        @test first(factors) == shared

        primitive = a * z + b
        monic(p) = divexact(p, leading_coefficient(p))
        @test Set(monic.(StructuralIdentifiability.fast_factor(shared * primitive))) ==
            Set(monic.([x, y, x + y, primitive]))
    end
end
