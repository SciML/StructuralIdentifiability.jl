@testset "Primality check (zerodim subroutine)" begin
    R, (x, y) = Nemo.PolynomialRing(Nemo.QQ, ["x", "y"])

    @test check_primality_zerodim([x^2 - 1, y^2 - 4]) == false

    @test check_primality_zerodim([(x + 5) * (x^3 - 7), y - 3]) == false

    @test check_primality_zerodim([x^3 - 5, y - 1]) == true

    @test check_primality_zerodim([x^2 + 1, y^3 - 3 * x + x + 5]) == true
end
