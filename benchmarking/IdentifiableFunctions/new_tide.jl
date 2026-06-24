using StructuralIdentifiability, StatsBase
using Logging

# macro taken from here
# https://discourse.julialang.org/t/help-writing-a-timeout-macro/16591/7
macro timeout(seconds, expr, fail)
    quote
        tsk = @task $(esc(expr))
        schedule(tsk)
        Timer($(esc(seconds))) do timer
            istaskdone(tsk) || Base.throwto(tsk, InterruptException())
        end
        try
            fetch(tsk)
        catch _
            $(esc(fail))
        end
    end
end

function simplify(genset)
    rff = StructuralIdentifiability.RationalFunctionField(genset)
    genset = StructuralIdentifiability.simplified_generating_set(rff, simplify=:weak)
    genset
end

function iterative_lie(
        ode::StructuralIdentifiability.ODE{P};
        initial_genset=nothing,
        tides=2,
        percentile=20,
        wait_for=20) where {P}
    genset = Vector()
    if initial_genset === nothing
        for y in ode.y_vars
            append!(genset, StructuralIdentifiability.extract_coefficients_ratfunc(ode.y_equations[y], ode.u_vars))
        end
    else
        genset = map(gen -> StructuralIdentifiability.parent_ring_change(gen, parent(ode)), initial_genset)
    end

    genset = simplify(genset)

    dense_id = frac -> length(numerator(frac)) * length(denominator(frac))
    dense_elems = empty(genset)
    new_elems = empty(genset)
    wave = 1
    while true
        @warn "Wave $wave" genset
        println(replace("$genset", "(t)" => ""))

        for gen in genset
            append!(new_elems, StructuralIdentifiability.extract_coefficients_ratfunc(
                StructuralIdentifiability.lie_derivative(gen, ode), ode.u_vars))
        end
        
        for tide in 1:tides
            @warn "The tide is coming: also computing derivatives of derivatives"
            new_new_elems = empty(new_elems)
            for gens in new_elems
                append!(new_new_elems, StructuralIdentifiability.extract_coefficients_ratfunc(
                    StructuralIdentifiability.lie_derivative(gens, ode), ode.u_vars))
            end
            append!(new_elems, new_new_elems)
        end
        filter!(!StructuralIdentifiability.is_rational_func_const, new_elems)
        
        @warn "Finding the simplest among $(length(new_elems)) new elems and $(length(dense_elems)) dense elems"
        new_elems = vcat(new_elems, dense_elems)
        sort!(new_elems, by=dense_id)
        containment = StructuralIdentifiability.field_contains(
            StructuralIdentifiability.RationalFunctionField(genset), new_elems, 0.99)
        @warn "Containment:" length(containment) count(iszero, containment)

        perm = filter(i -> containment[i] == 0, 1:length(containment))
        new_elems = new_elems[perm]
        n = length(new_elems)

        @warn "Densities:\n$(map(dense_id, new_elems))"
        if n == 0
            break
        end

        threshold = StatsBase.percentile(map(dense_id, new_elems), percentile)
        threshold_idx = findlast(map(dense_id, new_elems) .<= threshold)
        @warn """Simplifying elements 1:$threshold_idx ($percentile% percentile). 
                 Putting elems at $(threshold_idx+1)..$n back to dense list"""
        dense_elems = new_elems[threshold_idx+1:end]
        new_elems = new_elems[1:threshold_idx]

        new_genset = vcat(genset, new_elems)
        new_genset = @timeout wait_for begin 
            simplify(new_genset) 
        end Set()
        if isempty(new_genset)
            @warn "Interrupted by timeout"
            break
        end
        genset = new_genset
        wave += 1
    end

    return genset
end

TumorPillis2007 = StructuralIdentifiability.@ODEmodel(
            T'(t) =
                a * T(t) * (1 - b * T(t)) - c1 * N(t) * T(t) - D(t) * T(t) -
                KT * M(t) * T(t), #tumor cells
            L'(t) =
                -m * L(t) - q * L(t) * T(t) - ucte * L(t)^2 +
                r2 * C(t) * T(t) +
                pI * L(t) * I(t) / (gI + I(t)) +
                u1(t) - KL * M(t) * L(t), # tumor-specific effector cells, T-celss
            N'(t) =
                alpha1 - f * N(t) + g * T(t) * N(t) / (h + T(t)) - p * N(t) * T(t) -
                KN * M(t) * N(t), # non-specific effector cells, NK cells
            C'(t) = alpha2 - beta * C(t) - KC * M(t) * C(t), #circulating lymphocytes
            I'(t) =
                pt * T(t) * L(t) / (gt + T(t)) + w * L(t) * I(t) - muI * I(t) + u2(t), # IL-2, VI = u2 aplicación directa, terapia de IL2
            M'(t) = -gamma * M(t) + u1(t), #chemotherapy drug, terapia/aplicación de quimio, u1 = VM
            y1(t) = L(t),
            y2(t) = N(t),
            y3(t) = M(t)
        )
    

ode = TumorPillis2007

#find_ioequations(ode; loglevel=Logging.Debug)

genset = iterative_lie(ode, tides=2, percentile=20, wait_for=2)
println(genset)

