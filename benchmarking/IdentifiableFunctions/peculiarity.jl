using Nemo
using StructuralIdentifiability

R, (k01, k31, k21, k12, k13, k14, k41) =
    QQ[["k01", "k31", "k21", "k12", "k13", "k14", "k41"]...]

# Note: 7 elements
f1 = [
    k01,
    k12 + k13 + k14,
    k31 + k21 + k41,
    (k12^2 * k13 - k12^2 * k14 - k12 * k13^2 + k12 * k14^2 + k13^2 * k14 - k13 * k14^2) // (k31 * k12 - k31 * k14 - k21 * k13 + k21 * k14 - k12 * k41 + k13 * k41),
    (
        k31 * k12 * k13 - k31 * k13 * k14 - k21 * k12 * k13 + k21 * k12 * k14 -
        k12 * k14 * k41 + k13 * k14 * k41
    ) // (k31 * k12 - k31 * k14 - k21 * k13 + k21 * k14 - k12 * k41 + k13 * k41),
    (
        k31 * k21 * k12 - k31 * k21 * k13 + k31 * k13 * k41 - k31 * k14 * k41 -
        k21 * k12 * k41 + k21 * k14 * k41
    ) // (k31 * k12 - k31 * k14 - k21 * k13 + k21 * k14 - k12 * k41 + k13 * k41),
    (k31^2 * k21 - k31^2 * k41 - k31 * k21^2 + k31 * k41^2 + k21^2 * k41 - k21 * k41^2) // (k31 * k12 - k31 * k14 - k21 * k13 + k21 * k14 - k12 * k41 + k13 * k41),
]

# Note: 8 elements
f2 = [
    k01,
    k12 * k13 * k14,
    k31 * k21 * k41,
    k12 + k13 + k14,
    k31 + k21 + k41,
    k12 * k13 + k12 * k14 + k13 * k14,
    k31 * k13 + k21 * k12 + k14 * k41,
    k31 * k21 + k31 * k41 + k21 * k41,
]

# Prints: 1, 1
i12 = StructuralIdentifiability.check_field_membership(f1, f2, 0.9999)
i21 = StructuralIdentifiability.check_field_membership(f2, f1, 0.9999)
@info "" i12 i21

# Prints: 0, 0, ..., 0
res = fill(false, length(f2))
for i in 1:length(f2)
    all_except_i = f2[setdiff(1:length(f2), i)]
    res[i] = StructuralIdentifiability.check_field_membership(all_except_i, [f2[i]], 0.9999)
end
@info "" res
