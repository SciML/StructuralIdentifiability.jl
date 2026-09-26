using Test

@testset "Tests workflow reference" begin
    workflow = read(joinpath(@__DIR__, "..", ".github", "workflows", "Tests.yml"), String)
    @test occursin("grouped-tests.yml@v1", workflow)
    @test !occursin("grouped-tests.yml@1\"", workflow)
end
