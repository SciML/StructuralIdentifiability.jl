@testset "Constructive field membership" failfast = true begin
    R, (x,) = PolynomialRing(Nemo.QQ, ["x"])

    generators = [x^2, x^3]
    to_be_reduced = [x^2, x, 3one(R), zero(R)]

    memberships, remainders, relations_between_tags, tag_to_gen =
        StructuralIdentifiability.check_constructive_field_membership(
            StructuralIdentifiability.RationalFunctionField(generators),
            map(f -> f // one(f), to_be_reduced),
        )
    tags = gens(base_ring(parent(first(remainders))))

    @test length(tags) == 2
    @test all(memberships)
    @test map(string, remainders) == ["T1", "T1^2//T2", "3", "0"]
    @test tag_to_gen == Dict(tags[1] => x^2, tags[2] => x^3)
    @test length(relations_between_tags) == 1
    @test string(relations_between_tags[1]) == "T1^3 - T2^2"

    cases = []

    R, (T1,) = PolynomialRing(Nemo.QQ, ["T1"])
    append!(
        cases,
        [(generators = [T1^2], to_be_reduced = [T1, T1^2], memberships = Bool[0, 1])],
    )

    R, (x,) = PolynomialRing(Nemo.QQ, ["x"])
    append!(
        cases,
        [
            (
                generators = [(x - 1) // R(1), R(1) // (x^5 - 1), x // R(1)],
                to_be_reduced = [
                    (x^4 + x^3 + x^2 + x + 1) // one(R),
                    x // R(1),
                    R(33) // x^2,
                ],
                memberships = Bool[1, 1, 1],
            ),
            (
                generators = [(x^10 + x^9 + x^2 + 1) // (x^7 - x^6 - x^3 + 1)],
                to_be_reduced = [x // one(R), 2x // one(R), -3x // one(R)],
                memberships = Bool[0, 0, 0],
            ),
            (generators = [x^2], to_be_reduced = [x, x^88], memberships = Bool[0, 1]),
        ],
    )

    R, (x, y, z) = PolynomialRing(Nemo.QQ, ["x", "y", "z"])
    append!(
        cases,
        [
            (generators = [x, y], to_be_reduced = [x^2 + y^2, z], memberships = Bool[1, 0]),
            (
                generators = [x^2 + y^2, x^3 + y^3, x^4 + y^4],
                to_be_reduced = [x * y, x + y, x + y + 1, x + y + z],
                memberships = Bool[1, 1, 1, 0],
            ),
            (
                generators = [(x + y + z)^2, (x + y + z)^3, (x + y + z)^4],
                to_be_reduced = [(x + y + z)^18, x + 1, y + 2, z + 3],
                memberships = Bool[1, 0, 0, 0],
            ),
        ],
    )

    # NOTE: in this case it actually matter to cancel out the gcd after
    # computing the normal forms
    R, (a, b, y, x2, c, x1) = PolynomialRing(Nemo.QQ, ["a", "b", "y", "x2", "c", "x1"])
    append!(
        cases,
        [
            (
                generators = [
                    x1 // one(R),
                    a // one(R),
                    (a * c + c^2) // one(R),
                    c // x2,
                    x2 // (a + b),
                ],
                to_be_reduced = [
                    (a * c + c^2 + x1) // (a * c + c^2),
                    (a * c + c^2 + x1) // (a^2 + a * b + a * c + b * c),
                    (a * x2 + a * x1 + b * x1) // x2,
                ],
                memberships = Bool[1, 1, 1],
            ),
        ],
    )

    for case in cases
        generators = case.generators
        to_be_reduced = case.to_be_reduced
        memberships, remainders, relations_between_tags, tag_to_gen =
            StructuralIdentifiability.check_constructive_field_membership(
                StructuralIdentifiability.RationalFunctionField(generators),
                map(f -> f // one(f), to_be_reduced),
            )
        @test memberships == case.memberships
        tags = gens(base_ring(parent(first(remainders))))
        evaluate_tags = poly -> evaluate(poly, [tag_to_gen[tag] for tag in tags])
        for i in 1:length(relations_between_tags)
            @test iszero(evaluate_tags(relations_between_tags[i]))
        end
        for i in 1:length(remainders)
            if !memberships[i]
                continue
            end
            @test iszero(evaluate_tags(remainders[i]) - to_be_reduced[i])
        end
    end
end
