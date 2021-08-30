with(DifferentialAlgebra):
ring_diff := DifferentialRing(blocks = [[x15, x13, x7, x9, x2, x8, x14, x5, x4, x12, x6, x1, x10, x11, x3], [y4, y6, y2, y1, y3, y5] , [u]], derivations = [t]):
sys := [
diff(x15(t), t) - (-1/2500*x15(t) + 1/2000000*x7(t)),
diff(x13(t), t) - (e2a*x14(t) + 1/2*x10(t)*x6(t) - x2(t)*x13(t) - 1/50000*x13(t)),
diff(x7(t), t) - (5*i1*x6(t) - 1/2*x11(t)*x7(t)),
diff(x9(t), t) - (-1/2500*x9(t) - 1/2000000*x7(t)),
diff(x2(t), t) - (-1/5*x10(t)*x2(t) - x2(t)*x13(t) - x2(t)*k_deg - x2(t)*k2*u(t)*x8(t) - x2(t)*k3 + x1(t)*u(t)*k1 + x4(t)*t1 + x5(t)*t2),
diff(x8(t), t) - (1/2*x9(t) - c5*x8(t)),
diff(x14(t), t) - (-5*e2a*x14(t) + 1/2*x11(t)*x7(t)),
diff(x5(t), t) - (x2(t)*x13(t) - x5(t)*t2),
diff(x4(t), t) - (1/5*x10(t)*x2(t) - x4(t)*t1),
diff(x12(t), t) - (1/2000000*x7(t) - c3a*x12(t)),
diff(x6(t), t) - (-1/2*x10(t)*x6(t) + 1/50000*x13(t) - i1*x6(t) - x5(t)*t2),
diff(x1(t), t) - (-x1(t)*k_deg - x1(t)*u(t)*k1 + k_prod),
diff(x10(t), t) - (-1/5*x10(t)*x2(t) - x10(t)*i1a - 1/2*x10(t)*x6(t) - 1/10000*x10(t) + 1/2000*x11(t) + x12(t)*c4a),
diff(x11(t), t) - (5*x10(t)*i1a - 1/2*x11(t)*x7(t) - 1/400*x11(t)),
diff(x3(t), t) - (x2(t)*k2*u(t)*x8(t) + x2(t)*k3 - x3(t)*k_deg),
y4(t) - (x2(t) + x3(t) + x1(t)),
y6(t) - (x12(t)),
y2(t) - (x10(t) + x13(t)),
y1(t) - (x7(t)),
y3(t) - (x9(t)),
y5(t) - (x2(t))
];
res := CodeTools[CPUTime](RosenfeldGroebner(sys, ring_diff, singsol=none));