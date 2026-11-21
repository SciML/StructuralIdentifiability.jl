"""
    cmp_weighted(ode, W)

Comparator with arbitrary weights.
"""
function cmp_weighted(ode, W::Matrix{Int})
    @assert all(>=(0), W) "Weights must be non-negative."
    states, params = ode.x_vars, ode.parameters
    df(f) = lie_derivative(f, ode)
    features(f) = [
        !iszero(f),
        total_degree(denominator(f)),
        total_degree(numerator(f)),
        length(numerator(f)) + length(denominator(f)),
        !issubset(vars(f), params),
        total_degree_wrt(numerator(f), states),
        total_degree_wrt(denominator(f), states),
        number_of_terms(f, states),

        # Lie derivative
        !iszero(df(f)),
        total_degree(numerator(df(f))),
        total_degree(denominator(df(f))),
        length(numerator(df(f))) + length(denominator(df(f))),
        total_degree_wrt(df(f), states),
        total_degree_wrt(denominator(df(f)), states),
        number_of_terms(df(f), states),
    ]

    fallback = rational_function_cmp
    cmp = function (f, g)
        f = parent_ring_change(f, parent(ode))
        g = parent_ring_change(g, parent(ode))
        v_f = features(f)
        v_g = features(g)
        score_f = W * v_f
        score_g = W * v_g
        score_f != score_g && return score_f < score_g
        return fallback(f, g)
    end
    return cmp
end

# Matrices of weights for existing named comparators.
const matrix_comparators = let
    row(indices...) = [Int(j in indices) for _ in 1:1, j in 1:15]
    ordinary = vcat(row(1), row(2, 3), row(4), row(2))
    state = vcat(row(6, 7), row(7), row(8))
    lie = vcat(row(9), row(10, 11), row(12), row(11))
    lie_state = vcat(row(13), row(14), row(15))
    parameter_first = vcat(row(5), ordinary)

    matrices = Dict(
        cmp_default => parameter_first,
        cmp_lie => lie,
        cmp_states => vcat(state, ordinary),
        cmp_lie_states => vcat(lie_state, lie),
        cmp_pushka => vcat(state, lie_state, lie),
    )
    nrows = maximum(size(W, 1) for W in values(matrices))
    Dict(cmp => vcat(W, zeros(Int, nrows - size(W, 1), 15)) for (cmp, W) in matrices)
end
