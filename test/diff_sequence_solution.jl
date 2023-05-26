@testset "Computing variations around a sequence solution" begin
    function use_lie_derivatives(
        dds::ODE{P},
        params::Dict{P, T},
        ic::Dict{P, T},
        inputs::Dict{P, Array{T, 1}},
        num_terms::Int,
    ) where {T <: Generic.FieldElem, P <: MPolyElem{T}}
        explicit_sol = merge(
            Dict(x => Vector{Any}([x]) for (x, eq) in dds.x_equations),
            Dict(y => Vector{Any}([eq]) for (y, eq) in dds.y_equations)
        )
        for i in 2:num_terms
            for (k, v) in explicit_sol
                push!(explicit_sol[k], eval_at_dict(v[end], dds.x_equations))
            end
        end
        part_diffs = Dict((f, p) => [] for f in keys(explicit_sol) for p in vcat(dds.x_vars, dds.parameters))
        for i in 1:num_terms
            eval_dict = merge(
                params, ic,
                Dict(u => val[i] for (u, val) in inputs)
            )
            for (f, ders) in explicit_sol
                for p in vcat(dds.x_vars, dds.parameters)
                    push!(part_diffs[(f, p)], eval_at_dict(derivative(ders[i], p), eval_dict))
                end
            end
        end
        return part_diffs
    end

    locQQ = StructuralIdentifiability.Nemo.QQ

    ode = @ODEmodel(
        a'(t) = (23 * k1 * a(t) - 7 * b(t)^3) // (a(t)^2 + b(t)^2) - c(t)^3 * k1 * b(t),
        b'(t) = a(t) + 17 * (b(t) - c(t))^2 + 1 // (a(t) + b(t) - k2),
        y(t) = (a(t) + b(t) - c(t)) // (a(t)^2 + k2)
    )
    params = Dict(k1 => locQQ(1), k2 => locQQ(2))
    ic = Dict(a => locQQ(3), b => locQQ(-4))
    inputs = Dict(c => [locQQ(5), locQQ(-6), locQQ(7), locQQ(-8)])
    seq_sol, diff_sol = differentiate_sequence_solution(ode, params, ic, inputs, 2)
    diff_y = differentiate_sequence_output(ode, params, ic, inputs, 2)
    lie_ders_sol = use_lie_derivatives(ode, params, ic, inputs, 2)
    merged = merge(diff_sol, diff_y)
    for k in keys(lie_ders_sol)
        @info k
        @test merged[k] == lie_ders_sol[k]
    end
end
