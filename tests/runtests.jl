using Test
using TestSetExtensions

using Oscar

include("../io_equation.jl")
include("../power_series_utils.jl")

function random_ps(ps_ring, range = 1000)
    result = zero(ps_ring)
    t = gen(ps_ring)
    for i in 0:(max_precision(ps_ring) - 1)
        result += (rand(Int) % range) * t^i
    end
    return result
end

function random_ps_matrix(ps_ring, matrix_space)
    return map(e -> random_ps(ps_ring), zero(matrix_space))
end


@info "Testing started"

@testset "All the tests" begin
    @includetests ARGS
end

