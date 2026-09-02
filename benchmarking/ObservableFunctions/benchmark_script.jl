using AbstractAlgebra
using DataStructures
using StructuralIdentifiability
using StructuralIdentifiability: lie_derivatives_up_to
using RationalFunctionFields
using PrettyTables
using BenchmarkTools

# -----------------------------------------------------------------------------
# This script benchmarks StructuralIdentifiability.jl's algorithm for computing
# identifiable functions against a "naive" baseline that expands Lie derivatives
# up to a brute-force-determined order and then simplifies the resulting
# rational function field directly. It reports, per model, the prolongation
# orders/degree/size needed by each approach and the resulting runtimes.
# -----------------------------------------------------------------------------

"""
    prolongation_dict(ode, ords)

Build the `Dict` mapping each output equation `ode.y_equations[y]` to its
prolongation order `ords[i]`, as expected by `lie_derivatives_up_to`.
"""
function prolongation_dict(ode, ords)
    return Dict(ode.y_equations[y] => ords[i] for (i, y) in enumerate(ode.y_vars))
end

"""
    get_prolongation_orders(ode, with_identifiable_functions) -> Vector{Int}

For each output of `ode`, find the smallest prolongation order (number of Lie derivatives) needed to
get the full field of observable functions.

If `with_identifiable_functions`, the computations strats with the field of 
identifiable functions in parameter.

Returns a vector of orders, one per entry of `ode.y_vars`.
"""
function get_prolongation_orders(ode, with_identifiable_functions)
    orders = [0 for _ in ode.y_vars]
    original_gens = [values(ode.y_equations)...]
    if with_identifiable_functions
        original_gens = vcat(original_gens, find_identifiable_functions(ode))
    end
    rff = RationalFunctionField(original_gens)
    curr_trdeg = trdeg(rff)

    # Independently grow the order for each output until the transcendence
    # degree of the accumulated field saturates.
    for i in 1:length(ode.y_vars)
        while true
            orders[i] += 1
            original_gens = vcat(original_gens, lie_derivatives_up_to(ode, prolongation_dict(ode, orders)))
            new_trdeg = trdeg(RationalFunctionField(original_gens))
            if new_trdeg == curr_trdeg
                break
            else
                curr_trdeg = new_trdeg
            end
        end
    end
    return orders
end

"""
    trdeg(rff) -> Int

Compute the transcendence degree of a `RationalFunctionField`
"""
function trdeg(rff)
    RationalFunctionFields.update_trbasis_info!(rff, 0.99)
    return length(rff.trbasis)
end

function max_over_list(funcs, property)
    return maximum(f -> property(denominator(f)) + property(numerator(f)), funcs)
end

"""
    get_deg_size(ode, ords) -> (deg, size)

Compute the Lie-derivative prolongations of `ode` up to orders `ords`
(one per output, see `prolongation_dict`) and return the worst-case total
degree and worst-case size (`deg`, `size`) among them.
"""
function get_deg_size(ode, ords)
    prolongations = lie_derivatives_up_to(ode, prolongation_dict(ode, ords))
    deg = max_over_list(prolongations, total_degree)
    size = max_over_list(prolongations, length)
    return (deg, size)
end

"""
    get_orders_and_degrees(ode) -> Vector

Summarize the prolongation cost of `ode` for the benchmarking table. Computes:
- `ord_low`/`deg_low`/`size_low`: prolongation orders (without identifiable
  functions as extra generators) and the resulting degree/size;
- `ord_high`/`deg_high`/`size_high`: the same but with identifiable functions
  included as extra generators;
- `deg_final`/`size_final`: degree/size of the model's identifiable functions
  themselves (including states), as returned by the actual algorithm.

Returns a flat vector in that order, ready to be stacked into a table row.
"""
function get_orders_and_degrees(ode)
    ord_low = get_prolongation_orders(ode, false)
    ord_high = get_prolongation_orders(ode, true)

    id_funcs = find_identifiable_functions(ode, with_states = true)

    return [
        ord_low, get_deg_size(ode, ord_low)...,
        ord_high, get_deg_size(ode, ord_high)...,
        max_over_list(id_funcs, total_degree),
        max_over_list(id_funcs, length)
    ]
end

"""
    compute_direct(ode)

Direct approach for computing observable functions: expand Lie derivatives
up to the orders found by `get_prolongation_orders(ode, false)`, then
directly simplify the generating set of the resulting rational function field.
"""
function compute_direct(ode)
    ords = get_prolongation_orders(ode, false)
    prolongations = lie_derivatives_up_to(ode, prolongation_dict(ode, ords))
    return simplified_generating_set(RationalFunctionField(prolongations))
end

"""
    get_runtimes(models, to_skip) -> Vector

Benchmark, for each `(name, ode)` pair in `models`, the runtime of
`find_identifiable_functions` (the algorithm) versus `compute_direct`. 
Models listed in `to_skip` (by name) skip the naive benchmark
(e.g. because it is prohibitively slow) and are recorded with `Inf` time and
`Inf` speedup instead.

Returns a vector of rows `[name, time_algo, time_naive, speedup]`, where
`speedup = time_naive / time_algo`.
"""
function get_runtimes(models, to_skip)
    result = []
    for (name, ode) in models
        time_algo = @belapsed find_identifiable_functions($ode, with_states=true)
        time_naive = Inf
        if !(name in to_skip)
            time_naive = @belapsed compute_direct($ode)
        end
        push!(result, [name, time_algo, time_naive, time_naive / time_algo])
    end
    return result
end

###########################
# Example models
###########################

siwr = @ODEmodel(
    S'(t) = -beta_W*S(t)*W(t)-beta_I*S(t)*I(t),
    I'(t) = beta_W*S(t)*W(t)+beta_I*S(t)*I(t) - gamma*I(t),
    W'(t) = alpha*I(t)-zeta*W(t),
    y(t) = W(t)
)

MM = @ODEmodel(
    x1'(t) = V1 * x1(t) / (1 + K1 * x1(t) * (1 + I(t) / L1) + K2 * x2(t) * (1 + I(t) / L2)),
    x2'(t) = V2 * x2(t) / (1 + K1 * x1(t) * (1 + I(t) / L1) + K2 * x2(t) * (1 + I(t) / L2)),
    y(t) = x1(t) + x2(t)
)

sliqr = @ODEmodel(
    S'(t) = -b * In(t) * S(t) / N - u(t) * S(t) / N,
    L'(t) = b * In(t) * S(t) / N - a * L(t),
    In'(t) = a * L(t) - g * In(t) + s * Q(t),
    Q'(t) = (1 - e) * g * In(t) - s * Q(t),
    y(t) = In(t) / N
)

cancer1 = @ODEmodel(
    x'(t) = mu_m * (1 - q / Q(t)) * x(t) - (v(t) * R / (R + Q(t)) + delta * x(t)) * x(t),
    v'(t) = -d * v(t),
    Q'(t) = (g1 * u(t) + g2) * (Q_m - Q(t)) - mu_m * (Q(t) - q),
    P'(t) = b * Q(t) + s * x(t) * Q(t) - e * P(t),
    y(t) = P(t), y2(t) = Q(t)
)

cancer2 = @ODEmodel(
    x'(t) = mu_m * (1 - q / Q(t)) * x(t) - (v(t) * R / (R + Q(t)) + delta * x(t)) * x(t),
    v'(t) = -d * v(t),
    Q'(t) = (g1 * u(t) + g2) * (Q_m - Q(t)) - mu_m * (Q(t) - q),
    P'(t) = b * Q(t) + s * x(t) * Q(t) - e * P(t),
    y(t) = P(t), y2(t) = v(t)
)

eaihrd = @ODEmodel(
    A'(t) = a * E(t) - r1 * A(t),
    I'(t) = s * E(t) - (h + r2) * I(t),
    H'(t) = h * I(t) - (r3 + d) * H(t),
    R'(t) = r1 * A(t) + r2 * I(t) + r3 * H(t),
    D'(t) = d * H(t),
    E'(t) =
        (N - A(t) - I(t) - H(t) - R(t) - D(t) - E(t)) * (c1 * A(t) + c2 * I(t)) - (a + s) * E(t),
    y(t) = D(t) 
)

models = OrderedDict(
    :SIWR => siwr,
    :MM => MM,
    :SLIQR => sliqr,
    :cancer1 => cancer1,
    :cancer2 => cancer2,
    :EAIHRD => eaihrd
)

###########################
# Run benchmarks and report
###########################

# Table 1: prolongation orders and degree/size, with and without identifiable
# functions as extra generators, plus the degree/size of the final identifiable
# functions themselves.
degsize_data = stack([[name, get_orders_and_degrees(ode)...] for (name, ode) in models]; dims = 1)

pretty_table(degsize_data; column_labels = ["Name", "OrdLow", "DegLow", "SizeLow", "OrdHigh", "DegHigh", "SizeHigh", "DegFinal", "SizeFinal"])

# Table 2: runtime of the algorithm vs. the direct baseline (EAIHRD is skipped
# for the naive computation as it is too slow).
runtimes = get_runtimes(models, [:EAIHRD])
pretty_table(stack(runtimes; dims = 1); column_labels = ["Name", "Time algo", "Time naive", "Speedup"])
