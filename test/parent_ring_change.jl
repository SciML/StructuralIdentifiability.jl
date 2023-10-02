@testset "Parent ring change" begin
    for field in [Nemo.QQ, Nemo.GF(2^31 - 1)]
        R, (z, x, y) = PolynomialRing(QQ, ["z", "x", "y"], ordering = :degrevlex)
        R_, (y_, x_) = PolynomialRing(QQ, ["y", "x"], ordering = :lex)
        R__, (x__, t__, y__, z__) =
            PolynomialRing(QQ, ["x", "t", "y", "z"], ordering = :deglex)

        f = 2x + 3y + x^7 * y
        @test f == StructuralIdentifiability.parent_ring_change(f, R, matching = :byname)
        @test f == StructuralIdentifiability.parent_ring_change(f, R, matching = :byindex)

        f_ = StructuralIdentifiability.parent_ring_change(f, R_, matching = :byname)
        f__ = StructuralIdentifiability.parent_ring_change(f, R__, matching = :byname)
        @test f_ == 2x_ + 3y_ + x_^7 * y_
        @test f__ == 2x__ + 3y__ + x__^7 * y__
        @test f == StructuralIdentifiability.parent_ring_change(f_, R, matching = :byname)
        @test f == StructuralIdentifiability.parent_ring_change(f__, R, matching = :byname)

        @test_throws ArgumentError StructuralIdentifiability.parent_ring_change(
            x + z,
            R_,
            matching = :byname,
        )

        f = 3y
        f_ = StructuralIdentifiability.parent_ring_change(f, R__, matching = :byname)
        @test f_ == 3y__

        f__ = 2x__ + 5t__^3 + 3y__^2
        f = StructuralIdentifiability.parent_ring_change(f__, R, matching = :byindex)
        @test_throws ArgumentError StructuralIdentifiability.parent_ring_change(
            f__,
            R_,
            matching = :byindex,
        )
        @test f == 2z + 5x^3 + 3y^2
        @test f__ ==
              StructuralIdentifiability.parent_ring_change(f, R__, matching = :byindex)

        R, (x,) = PolynomialRing(field, ["x"])
        R2, (x2,) = PolynomialRing(Nemo.FractionField(R), ["x"])
        R3, (x3,) = PolynomialRing(field, ["x"])
        f = x3^2
        f_ = StructuralIdentifiability.parent_ring_change(f, R2)
        @test parent(f) == R3
        @test parent(f_) == R2
        @test repr(f_) == "x^2"
    end
end
