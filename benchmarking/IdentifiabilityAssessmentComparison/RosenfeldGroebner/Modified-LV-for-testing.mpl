with(DifferentialAlgebra):
ring_diff := DifferentialRing(blocks = [[x1, x2], [y1] ], derivations = [t]):
sys := [
diff(x1(t), t) - (a*x1(t) + b*x1(t) - x2(t)*c*x1(t)),
diff(x2(t), t) - (-a*b*x2(t) + d*x2(t)*x1(t)),
y1(t) - (x1(t))
];
res := CodeTools[CPUTime](RosenfeldGroebner(sys, ring_diff, singsol=none));