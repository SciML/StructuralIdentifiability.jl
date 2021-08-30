with(DifferentialAlgebra):
ring_diff := DifferentialRing(blocks = [[x5, x6, x4, x2, x1, x3], [y1, y2] ], derivations = [t]):
sys := [
diff(x5(t), t) - (k5*x6(t) + k4*x6(t) - k6*x5(t)*x3(t)),
diff(x6(t), t) - (-k5*x6(t) - k4*x6(t) + k6*x5(t)*x3(t)),
diff(x4(t), t) - (-k3*x4(t) - k2*x4(t) + k1*x1(t)*x2(t)),
diff(x2(t), t) - (k3*x4(t) + k2*x4(t) + k1*x1(t)*x2(t)),
diff(x1(t), t) - (k4*x6(t) + k2*x4(t) - k1*x1(t)*x2(t)),
diff(x3(t), t) - (k5*x6(t) + k3*x4(t) - k6*x5(t)*x3(t)),
y1(t) - (x3(t)),
y2(t) - (x2(t))
];
res := CodeTools[CPUTime](RosenfeldGroebner(sys, ring_diff, singsol=none));