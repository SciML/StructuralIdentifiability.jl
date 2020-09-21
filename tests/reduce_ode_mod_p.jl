@testset "Reducing ODE mod p" begin
    
    R, (x, y) = PolynomialRing(QQ, ["x", "y"])
    T = typeof(x)
    ode = ODE{T}(
        Dict{T, Union{T, Generic.Frac{T}}}(x => 1//2 * x - 5 * y, y => (3//5 * x + 17 * y) // (y - 1)), 
        Array{T, 1}()
    )

    ode_red = reduce_ode_mod_p(ode, 17)

    x, y = gens(ode_red.poly_ring)
    @test ode_red.equations[x] == 9 * x + 12 * y
    @test ode_red.equations[y] == 4 * x // (y + 16)
    
end
