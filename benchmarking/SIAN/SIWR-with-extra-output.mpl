read "../IdentifiabilityODE.mpl";

sys := [
diff(S(t), t) = -bi*S(t)*I(t) - S(t)*mu - S(t)*bw*W(t) + R(t)*a + mu,
diff(I(t), t) = bi*S(t)*I(t) - gam*I(t) + S(t)*bw*W(t) - mu*I(t),
diff(W(t), t) = xi*I(t) - xi*W(t),
diff(R(t), t) = gam*I(t) - R(t)*mu - R(t)*a,
y(t) = k*I(t),
y2(t) = S(t) + R(t) + I(t)
];
CodeTools[CPUTime](IdentifiabilityODE(sys, GetParameters(sys)));