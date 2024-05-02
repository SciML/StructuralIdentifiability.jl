using StructuralIdentifiability, StatsBase
using Logging
using Nemo


global_logger(SimpleLogger(stdout, Logging.Debug))

StructuralIdentifiability._groebner_loglevel[] = 10

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
        @warn "Checking containment"
        @warn "Gens: $genset"
        @warn "To check $new_elems"
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
        @warn "Simplifying $new_genset"
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
    

MAPK = @ODEmodel(
    KS00'(t) =
        -a00 * K(t) * S00(t) +
        b00 * KS00(t) +
        gamma0100 * FS01(t) +
        gamma1000 * FS10(t) +
        gamma1100 * FS11(t),
    KS01'(t) =
        -a01 * K(t) * S01(t) + b01 * KS01(t) + c0001 * KS00(t) - alpha01 * F(t) * S01(t) +
        beta01 * FS01(t) +
        gamma1101 * FS11(t),
    KS10'(t) =
        -a10 * K(t) * S10(t) + b10 * KS10(t) + c0010 * KS00(t) - alpha10 * F(t) * S10(t) +
        beta10 * FS10(t) +
        gamma1110 * FS11(t),
    FS01'(t) =
        -alpha11 * F(t) * S11(t) +
        beta11 * FS11(t) +
        c0111 * KS01(t) +
        c1011 * KS10(t) +
        c0011 * KS00(t),
    FS10'(t) = a00 * K(t) * S00(t) - (b00 + c0001 + c0010 + c0011) * KS00(t),
    FS11'(t) = a01 * K(t) * S01(t) - (b01 + c0111) * KS01(t),
    K'(t) = a10 * K(t) * S10(t) - (b10 + c1011) * KS10(t),
    F'(t) = alpha01 * F(t) * S01(t) - (beta01 + gamma0100) * FS01(t),
    S00'(t) = alpha10 * F(t) * S10(t) - (beta10 + gamma1000) * FS10(t),
    S01'(t) =
        alpha11 * F(t) * S11(t) - (beta11 + gamma1101 + gamma1110 + gamma1100) * FS11(t),
    S10'(t) =
        -a00 * K(t) * S00(t) + (b00 + c0001 + c0010 + c0011) * KS00(t) -
        a01 * K(t) * S01(t) + (b01 + c0111) * KS01(t) - a10 * K(t) * S10(t) +
        (b10 + c1011) * KS10(t),
    S11'(t) =
        -alpha01 * F(t) * S01(t) + (beta01 + gamma0100) * FS01(t) -
        alpha10 * F(t) * S10(t) + (beta10 + gamma1000) * FS10(t) -
        alpha11 * F(t) * S11(t) +
        (beta11 + gamma1101 + gamma1110 + gamma1100) * FS11(t),
    y1(t) = F(t),
    y2(t) = S00(t),
    y3(t) = S01(t) + S10(t),
    y4(t) = S11(t)
)

NFkB = @ODEmodel(
    x1'(t) = k_prod - k_deg * x1(t) - k1 * x1(t) * u(t),
    x2'(t) =
        -k3 * x2(t) - k_deg * x2(t) - a2 * x2(t) * x10(t) + t1 * x4(t) -
        a3 * x2(t) * x13(t) +
        t2 * x5(t) +
        (k1 * x1(t) - k2 * x2(t) * x8(t)) * u(t),
    x3'(t) = k3 * x2(t) - k_deg * x3(t) + k2 * x2(t) * x8(t) * u(t),
    x4'(t) = a2 * x2(t) * x10(t) - t1 * x4(t),
    x5'(t) = a3 * x2(t) * x13(t) - t2 * x5(t),
    x6'(t) = c6a * x13(t) - a1 * x6(t) * x10(t) + t2 * x5(t) - i1 * x6(t),
    x7'(t) = i1 * kv * x6(t) - a1 * x11(t) * x7(t),
    x8'(t) = c4 * x9(t) - c5 * x8(t),
    x9'(t) = c2 - c1 * x7(t) - c3 * x9(t),
    x10'(t) =
        -a2 * x2(t) * x10(t) - a1 * x10(t) * x6(t) + c4a * x12(t) - c5a * x10(t) -
        i1a * x10(t) + e1a * x11(t),
    x11'(t) = -a1 * x11(t) * x7(t) + i1a * kv * x10(t) - e1a * kv * x11(t),
    x12'(t) = c2a + c1a * x7(t) - c3a * x12(t),
    x13'(t) = a1 * x10(t) * x6(t) - c6a * x13(t) - a3 * x2(t) * x13(t) + e2a * x14(t),
    x14'(t) = a1 * x11(t) * x7(t) - e2a * kv * x14(t),
    x15'(t) = c2c + c1c * x7(t) - c3c * x15(t),
    y1(t) = x7(t),
    y2(t) = x10(t) + x13(t),
    y3(t) = x9(t),
    y4(t) = x1(t) + x2(t) + x3(t),
    y5(t) = x2(t),
    y6(t) = x12(t),
)

# surprisingly hard for the tides but fine for IO
JAKS = @ODEmodel(
	x1'(t) = t6*x2(t) - t5*x1(t) - 2*u(t)*x1(t)*t1,
	x2'(t) = -t6*x2(t) + t5*x1(t),
	x3'(t) = x6(t)*x3(t)*t2 - 3*x3(t)*t2 + 2*u(t)*x1(t)*t1,
	x4'(t) = -t3*x4(t) - x6(t)*x3(t)*t2 + 3*x3(t)*t2,
	x5'(t) = t3*x4(t) - x5(t)*t4,
	x6'(t) = (-x6(t)*x3(t)*x10(t)*t7*t13 - x6(t)*x3(t)*t7 - 92*x6(t)*x10(t)*x1(t)*t8*t13^2 - 92*x6(t)*x10(t)*t8*t13 - 92*x6(t)*x1(t)*t8*t13 - x6(t)*x1(t)*t7*t13*x4(t) - 92*x6(t)*t8 - x6(t)*t7*x4(t) + 276*x10(t)*x1(t)*t8*t13^2 + 276*x10(t)*t8*t13 + 276*x1(t)*t8*t13 + 276*t8)//(x10(t)*x1(t)*t13^2 + x10(t)*t13 + x1(t)*t13 + 1),
	x7'(t) = -92*x7(t)*t10 + x7(t)*x6(t)*t9 - 3*x7(t)*t9 + 15180*t10,
	x8'(t) = -x7(t)*t11 + 165*t11,
	x9'(t) = -2*x9(t)*u(t)*t12,
	x10'(t) = (-x8(t)*t16*x10(t) + x8(t)*t14 - t16*x10(t)*t15)//(x8(t) + t15),
	y1(t) = x3(t) + x1(t) + x4(t),
	y2(t) = -x9(t)*t18 + x5(t)*t18 + t18*x3(t) + t18*x4(t) + 1//3*t18,
	y3(t) = t19*x5(t) + t19*x4(t),
	y4(t) = -t20*x6(t) + 3*t20,
	y5(t) = x8(t)*t21,
	y6(t) = (x8(t)*t22*t17)//t11,
	y7(t) = x10(t),
	y8(t) = -x7(t) + 165
)

ode = NFkB


@warn collect(gens(ode.poly_ring))

#find_identifiable_functions(ode; with_states=true, loglevel=Logging.Debug)

genset = iterative_lie(ode, tides=1, percentile=20, wait_for=3000)
println(genset)

