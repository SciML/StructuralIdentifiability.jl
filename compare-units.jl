using Revise, StructuralIdentifiability, ParamPunPam, RationalFunctionFields, Unitful, Nemo, AbstractAlgebra
import StructuralIdentifiability: parent_ring_change, lie_derivative
import RationalFunctionFields: fields_equal

function eval_at_dict(poly::P, d::Dict{P,S}) where {P<:MPolyRingElem,S}
    R = parent(poly)
    @assert R == parent(first(keys(d)))
    xs = gens(parent(first(keys(d))))
    xs_sym = [get(d, x, 0) for x in xs if string(x) in map(string, gens(R))]
    accum = nothing
    try
    for t in terms(poly)
        cf = coeff(t, 1)
        cf_ = BigInt(numerator(cf)) // BigInt(denominator(cf))
        ex = exponent_vector(t, 1)
        if accum == nothing
            accum = cf_ * prod(xs_sym .^ ex)
        else
            accum += cf_ * prod(xs_sym .^ ex)
        end
    end
    catch
        return "beda!"
    end 
    return accum
end

eval_at_dict(frac::AbstractAlgebra.Generic.FracFieldElem{T}, d) where {T} = 
   begin
    num = eval_at_dict(numerator(frac), d)
    den = eval_at_dict(denominator(frac), d)
    "beda!" in (num, den) && return "beda!"
    num // den
   end

ode = @ODEmodel(
    S'(t) = -beta_W*S(t)*W(t) - beta_I*S(t)*I(t),
    I'(t) = beta_W*S(t)*W(t) + beta_I*S(t)*I(t) - gamma*I(t),
    W'(t) = alpha*I(t) - zeta*W(t),
    y(t) = W(t)
)

# Domain knowledge is that the units of S and I are the same.
units = Dict(
    S => u"P",        # meaning: "population"
    I => u"P",       
    W => u"C",        # meaning: "concentration"
    gamma => u"1/s",  # meaning: "1/second"
    zeta => u"1/s",
    beta_W => u"1/C * 1/s", 
    beta_I => u"1/P * 1/s",
    alpha => u"C * 1/P * 1/s",
)

cmp = (f, g) -> begin
    f = parent_ring_change(f, parent(ode))
    g = parent_ring_change(g, parent(ode))
    uf = Unitful.unit(eval_at_dict(f, units))
    ug = Unitful.unit(eval_at_dict(g, units))
    length(string(uf)) < length(string(ug))
end

id1 = find_identifiable_functions(ode, with_states=true, cmp = cmp)

id2 = find_identifiable_functions(ode, with_states=true, cmp = RationalFunctionFields.rational_function_cmp)

@assert fields_equal(RationalFunctionField(id1), RationalFunctionField(id2), 0.99)

u1 = Unitful.unit.(map(f -> eval_at_dict(f, units), id1))
u2 = Unitful.unit.(map(f -> eval_at_dict(f, units), id2))
