@testset "Power Series Differentiation" begin
    
    T, t = PowerSeriesRing(QQ, 300, "t"; model=:capped_absolute)
    
    @test ps_diff(zero(T)) == zero(T)
    
    @test ps_diff(one(T)) == zero(T)
    
    @test ps_diff(gen(T)) == one(T)
    
    @test ps_diff(t^3 - t^4 * (1 // 8) + 5 * t^10) == 3 * t^2 - t^3 * (1 // 2) + 50 * t^9
    
    @test ps_diff( (1 + t)^1000 ) == truncate(1000 * (1 + t)^999, 299)
    
end
