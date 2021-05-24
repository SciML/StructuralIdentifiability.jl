read "../IdentifiabilityODE.mpl";

sys := [
diff(S(t), t) = -b*S(t)*Ninv(t)*A(t)*q - b*S(t)*Ninv(t)*I(t) - b*S(t)*Ninv(t)*J(t),
diff(A(t), t) = -E(t)*k*r + E(t)*k - A(t)*g1,
diff(E(t), t) = b*S(t)*Ninv(t)*A(t)*q + b*S(t)*Ninv(t)*I(t) + b*S(t)*Ninv(t)*J(t) - E(t)*k,
diff(C(t), t) = alpha*I(t),
diff(Ninv(t), t) = 0,
diff(J(t), t) = alpha*I(t) - g2*J(t),
diff(I(t), t) = -alpha*I(t) + E(t)*k*r - g1*I(t),
y2(t) = Ninv(t),
y(t) = C(t)
];
CodeTools[CPUTime](IdentifiabilityODE(sys, GetParameters(sys)));