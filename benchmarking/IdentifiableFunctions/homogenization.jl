using StructuralIdentifiability, Groebner, Nemo, ParamPunPam

Bilirubin2_io = @ODEmodel(
    x1'(t) =
        -(k21 + k31 + k41 + k01) * x1(t) + k12 * x2(t) + k13 * x3(t) + k14 * x4(t) + u(t),
    x2'(t) = k21 * x1(t) - k12 * x2(t),
    x3'(t) = k31 * x1(t) - k13 * x3(t),
    x4'(t) = k41 * x1(t) - k14 * x4(t),
    y1(t) = x1(t)
)

#######
# Dancing for rain in order to extract some form of the MQS ideal generators

funcs = find_identifiable_functions(Bilirubin2_io, with_states = true, strategy = (:gb,))

rff = StructuralIdentifiability.RationalFunctionField(funcs)
cfs = StructuralIdentifiability.beautiful_generators(rff)
rff = StructuralIdentifiability.RationalFunctionField(cfs)

K = GF(2^31 - 1)
mqs = rff.mqs
xs = Nemo.gens(Nemo.parent(mqs))
ParamPunPam.reduce_mod_p!(mqs, K)
point = [rand(K) for _ in 1:(length(xs) - 1)]
ideal_spec = StructuralIdentifiability.specialize_mod_p(mqs, point)

#######

# 5 ms
@time gb = groebner(ideal_spec, ordering = Groebner.DegRevLex());

# There is an existent possibility that this would not finish in two and a half lifetimes
# @time gb = groebner(ideal_spec, ordering = Groebner.Lex(), loglevel = -3);

hom_ideal_spec = StructuralIdentifiability.homogenize(ideal_spec);
# 100 ms
@time Groebner.groebner(hom_ideal_spec, ordering = Groebner.Lex());
