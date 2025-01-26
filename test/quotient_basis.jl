@testset "Quotient basis" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"], internal_ordering = :degrevlex)

    @test Set(quotient_basis([x, y])) == Set([R(1)])
    @test Set(quotient_basis([x + y, y])) == Set([R(1)])
    @test Set(quotient_basis([-2x + 8, 10y + 1])) == Set([R(1)])

    @test Set(quotient_basis([x^2 + 1, y^2 - x - 1])) == Set([R(1), y, x, x * y])
    @test Set(quotient_basis([8x^5 + 1, y - 1])) == Set([R(1), x, x^2, x^3, x^4])
    @test Set(quotient_basis([x^2, y^2, x * y])) == Set([R(1), y, x])
    @test Set(quotient_basis([x^3, y^3, x * y])) == Set([R(1), y, x, y^2, x^2])
    @test length(quotient_basis([8x^77 + 1, y^34 - 1])) == 77 * 34

    R, (x, y) = polynomial_ring(QQ, ["x", "y"], internal_ordering = :lex)
    @test Set(quotient_basis([x, y])) == Set([R(1)])
    @test Set(quotient_basis([x + y, y])) == Set([R(1)])
    @test Set(quotient_basis([-2x + 8, 10y + 1])) == Set([R(1)])

    @test Set(quotient_basis([x^2 + 1, y^2 - 1])) == Set([R(1), y, x, x * y])
    @test Set(quotient_basis([8x^5 + 1, y - 1])) == Set([R(1), x, x^2, x^3, x^4])
    @test Set(quotient_basis([x^2, y^2, x * y])) == Set([R(1), y, x])
    @test Set(quotient_basis([x^3, y^3, x * y])) == Set([R(1), y, y^2, x, x^2])
    @test length(quotient_basis([8x^77 + 1, y^34 - 1])) == 77 * 34

    R, (x, y) = polynomial_ring(QQ, ["x", "y"], internal_ordering = :deglex)

    @test Set(quotient_basis([x, y])) == Set([R(1)])
    @test Set(quotient_basis([x + y, y])) == Set([R(1)])
    @test Set(quotient_basis([-2x + 8, 10y + 1])) == Set([R(1)])

    @test Set(quotient_basis([x^2 + 1, y^2 - x - 1])) == Set([R(1), y, x, x * y])
    @test Set(quotient_basis([8x^5 + 1, y - 1])) == Set([R(1), x, x^2, x^3, x^4])
    @test Set(quotient_basis([x^2, y^2, x * y])) == Set([R(1), y, x])
    @test Set(quotient_basis([x^3, y^3, x * y])) == Set([R(1), y, x, y^2, x^2])
    @test length(quotient_basis([8x^77 + 1, y^34 - 1])) == 77 * 34

    @test_throws DomainError("Input does not define zerodimensional ideal") quotient_basis([
        x^2 + y - 3,
    ])
end
