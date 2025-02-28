@testset "Identifiability of linear compartment models" begin
    test_cases = []

    push!(
        test_cases,
        Dict(
            :graph => [Array{Int, 1}(), [1], [1]],
            :outputs => [1],
            :leaks => Array{Int, 1}(),
            :inputs => [1],
            :result => Dict((3, 1) => :locally, (2, 1) => :locally),
        ),
    )

    #--------------------------------------------------------------------------

    push!(
        test_cases,
        Dict(
            :graph => [Array{Int, 1}(), [1], [1, 2]],
            :outputs => [1],
            :leaks => Array{Int, 1}(),
            :inputs => [1],
            :result => Dict(
                (3, 1) => :nonidentifiable,
                (3, 2) => :nonidentifiable,
                (2, 1) => :locally,
            ),
        ),
    )

    #--------------------------------------------------------------------------

    push!(
        test_cases,
        Dict(
            :graph => [[2], [1, 3], [1, 2]],
            :outputs => [1],
            :leaks => Array{Int, 1}(),
            :inputs => [2],
            :result => Dict(
                (2, 3) => :nonidentifiable,
                (3, 1) => :nonidentifiable,
                (1, 2) => :locally,
                (3, 2) => :nonidentifiable,
                (2, 1) => :globally,
            ),
        ),
    )

    #--------------------------------------------------------------------------

    push!(
        test_cases,
        Dict(
            :graph => [[2], [1, 3], [1, 2]],
            :outputs => [2],
            :leaks => Array{Int, 1}(),
            :inputs => [2],
            :result => Dict(
                (2, 3) => :nonidentifiable,
                (3, 1) => :nonidentifiable,
                (1, 2) => :locally,
                (3, 2) => :nonidentifiable,
                (2, 1) => :nonidentifiable,
            ),
        ),
    )

    #--------------------------------------------------------------------------

    push!(
        test_cases,
        Dict(
            :graph => [[2, 3], [1, 3], [1, 2]],
            :outputs => [3],
            :leaks => [1],
            :inputs => [1, 2],
            :result => Dict(
                (1, 2) => :globally,
                (1, 3) => :globally,
                (2, 1) => :globally,
                (2, 3) => :globally,
                (3, 1) => :globally,
                (3, 2) => :globally,
                (1, 0) => :globally,
            ),
        ),
    )

    #--------------------------------------------------------------------------

    push!(
        test_cases,
        Dict(
            :graph => [[2, 3], [1, 3], [1, 2]],
            :outputs => [3],
            :leaks => Array{Int, 1}(),
            :inputs => [1, 2],
            :result => Dict(
                (1, 2) => :nonidentifiable,
                (1, 3) => :globally,
                (2, 1) => :nonidentifiable,
                (2, 3) => :globally,
                (3, 1) => :nonidentifiable,
                (3, 2) => :nonidentifiable,
            ),
        ),
    )

    #--------------------------------------------------------------------------

    push!(
        test_cases,
        Dict(
            :graph => [[2], [1, 3], [1, 2]],
            :outputs => [3],
            :leaks => [1],
            :inputs => [1, 3],
            :result => Dict(
                (1, 2) => :globally,
                (2, 1) => :globally,
                (2, 3) => :globally,
                (3, 1) => :globally,
                (3, 2) => :globally,
                (1, 0) => :globally,
            ),
        ),
    )

    #--------------------------------------------------------------------------

    push!(
        test_cases,
        Dict(
            :graph => [[2], [1, 3], [1, 2]],
            :outputs => [3],
            :leaks => Array{Int, 1}(),
            :inputs => [1, 3],
            :result => Dict(
                (1, 2) => :nonidentifiable,
                (2, 1) => :nonidentifiable,
                (2, 3) => :nonidentifiable,
                (3, 1) => :nonidentifiable,
                (3, 2) => :nonidentifiable,
            ),
        ),
    )

    #--------------------------------------------------------------------------

    for case in test_cases
        ode = linear_compartment_model(
            case[:graph],
            inputs = case[:inputs],
            outputs = case[:outputs],
            leaks = case[:leaks],
        )
        bring = ode.poly_ring
        correct = Dict{QQMPolyRingElem, Symbol}()
        for (e, id) in case[:result]
            correct[str_to_var("a_$(e[2])_$(e[1])", bring)] = id
        end
        result = assess_identifiability(ode, funcs_to_check = collect(keys(correct)))
        @test correct == result
    end
end
