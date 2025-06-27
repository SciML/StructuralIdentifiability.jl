with(DifferentialThomas):
with(Tools):
Ranking([t], [[S, I, W, R], [y] ]):
sys := [
diff(S(t), t) - (-bi*S(t)*II(t) - S(t)*mu - S(t)*bw*W(t) + R(t)*a + mu),
diff(II(t), t) - (bi*S(t)*II(t) - gam*II(t) + S(t)*bw*W(t) - mu*II(t)),
diff(W(t), t) - (xi*II(t) - xi*W(t)),
diff(R(t), t) - (gam*II(t) - R(t)*mu - R(t)*a),
y(t) - (k*II(t))
];
res := CodeTools[CPUTime](ThomasDecomposition(sys));