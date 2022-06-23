@testset "Trie for exponents vectors" begin

for _ in 1:10
    d = rand([i + 10 for i in 1:10])
    t = ExpVectTrie(d)
    push!(t, [0 for _ in 1:d])
    vects = [rand([i for i in 1:10], d) for _ in 1:20]
    for v in vects
        push!(t, v)
    end
    for v in vects
        diff, best = get_max_below(t, v)
        @test diff == 0
        @test best == v
    end

    svects = Set(vects)
    push!(svects, [0 for _ in 1:d])

    for _ in 1:20
        v = rand([i for i in 1:10], d)
        diff, best = get_max_below(t, v)
        if v in svects
            @test diff == 0
            @test v == best
        else
            @test (best in svects)
            @test diff == sum(v .- best)
        end
    end
end end
