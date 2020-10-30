@testset "Primality check (univariate subroutine)" begin
    
    R, x = PolynomialRing(QQ, "x")

    @test check_primality_univariate([x^2 - 1, x^2 - 4]) == false

    @test check_primality_univariate([(x + 5) * (x^3 - 7)]) == false

    @test check_primality_univariate([x^3 - 5]) == true

    @test check_primality_univariate([x^2 + 1, x^3 - 3 * x + x + 5]) == true

end
