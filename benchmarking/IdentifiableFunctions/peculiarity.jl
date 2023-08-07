using Nemo
using StructuralIdentifiability

R, (k01, k31, k21, k12, k13, k14, k41) =
    QQ[["k01", "k31", "k21", "k12", "k13", "k14", "k41"]...]

# Note: 7 elements
f1_lists = [
    [R(1), k01],
    [R(1), k12 + k13 + k14],
    [R(1), k31 + k21 + k41],
    [
        k31 * k12 - k31 * k14 - k21 * k13 + k21 * k14 - k12 * k41 + k13 * k41,
        k12^2 * k13 - k12^2 * k14 - k12 * k13^2 + k12 * k14^2 + k13^2 * k14 - k13 * k14^2,
    ],
    [
        k31 * k12 - k31 * k14 - k21 * k13 + k21 * k14 - k12 * k41 + k13 * k41,
        k31 * k12 * k13 - k31 * k13 * k14 - k21 * k12 * k13 + k21 * k12 * k14 -
        k12 * k14 * k41 + k13 * k14 * k41,
    ],
    [
        k31 * k12 - k31 * k14 - k21 * k13 + k21 * k14 - k12 * k41 + k13 * k41,
        k31 * k21 * k12 - k31 * k21 * k13 + k31 * k13 * k41 - k31 * k14 * k41 -
        k21 * k12 * k41 + k21 * k14 * k41,
    ],
    [
        k31 * k12 - k31 * k14 - k21 * k13 + k21 * k14 - k12 * k41 + k13 * k41,
        k31^2 * k21 - k31^2 * k41 - k31 * k21^2 + k31 * k41^2 + k21^2 * k41 - k21 * k41^2,
    ],
]
f1_fracs = [
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
f2_lists = [
    [R(1), k01],
    [R(1), k12 * k13 * k14],
    [R(1), k31 * k21 * k41],
    [R(1), k12 + k13 + k14],
    [R(1), k31 + k21 + k41],
    [R(1), k12 * k13 + k12 * k14 + k13 * k14],
    [R(1), k31 * k13 + k21 * k12 + k14 * k41],
    [R(1), k31 * k21 + k31 * k41 + k21 * k41],
]
f2_fracs = [
    (k01) // R(1),
    (k12 * k13 * k14) // R(1),
    (k31 * k21 * k41) // R(1),
    (k12 + k13 + k14) // R(1),
    (k31 + k21 + k41) // R(1),
    (k12 * k13 + k12 * k14 + k13 * k14) // R(1),
    (k31 * k13 + k21 * k12 + k14 * k41) // R(1),
    (k31 * k21 + k31 * k41 + k21 * k41) // R(1),
]

@assert length(f1_lists) == length(f1_fracs) == 7
@assert length(f2_lists) == length(f2_fracs) == 8

# Prints: true; true
p = 0.99999
i12 = StructuralIdentifiability.check_field_membership(f1_lists, f2_fracs, p)
i21 = StructuralIdentifiability.check_field_membership(f2_lists, f1_fracs, p)
@info "" i12 i21

# Prints: false, false, ..., false
res = [false for _ in 1:length(f2_lists)]
for i in 1:length(f2_lists)
    all_except_i = f2_lists[setdiff(1:length(f2_lists), i)]
    res[i] =
        StructuralIdentifiability.check_field_membership(all_except_i, [f2_fracs[i]], p)[1]
end
@info "" res
