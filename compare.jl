using Revise, StructuralIdentifiability, ParamPunPam, RationalFunctionFields
import StructuralIdentifiability: parent_ring_change, lie_derivative
import RationalFunctionFields: fields_equal

cmp = (f, g) -> begin
    f = parent_ring_change(f, parent(ode))
    g = parent_ring_change(g, parent(ode))
    df = lie_derivative(f, ode)
    dg = lie_derivative(g, ode)
    default_cmp = RationalFunctionFields.rational_function_cmp
    if iszero(df) && iszero(dg)
        default_cmp(f, g)
    else
        default_cmp(df, dg)
    end    
end

ode = @ODEmodel(
    S'(t) = -beta_W*S(t)*W(t)-beta_I*S(t)*I(t),
    I'(t) = beta_W*S(t)*W(t)+beta_I*S(t)*I(t) - gamma*I(t),
    W'(t) = alpha*I(t)-zeta*W(t),
    y(t) = W(t)
)

# ode = @ODEmodel(
#     x1'(t) = x1 + u(t),
#     y(t) = x1
# )

id1 = find_identifiable_functions(ode, with_states=true, cmp = cmp)

id2 = find_identifiable_functions(ode, with_states=true, cmp = RationalFunctionFields.rational_function_cmp)

@assert fields_equal(RationalFunctionField(id1), RationalFunctionField(id2), 0.99)
