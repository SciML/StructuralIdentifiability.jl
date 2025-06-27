with(DifferentialThomas):
with(Tools):
Ranking([t], [[x2, x4, x1, x3], [y] ]):
sys := [
diff(x2(t), t) - (alpha*x1(t) - beta*x2(t)),
diff(x4(t), t) - ((gama*sigma*x2(t)*x4(t) - delta*sigma*x3(t)*x4(t)) / (x3(t))),
diff(x1(t), t) - ((-b*c*x1(t) - b*x1(t)*x4(t) + 1) / (c + x4(t))),
diff(x3(t), t) - (gama*x2(t) - delta*x3(t)),
y(t) - (x1(t))
];
res := CodeTools[CPUTime](ThomasDecomposition(sys));