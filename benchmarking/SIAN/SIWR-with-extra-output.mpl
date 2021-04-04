k*I(t)diff(R(t), t) = gam*I(t) - R(t)*mu - R(t)*a,
diff(I(t), t) = bi*S(t)*I(t) - gam*I(t) + S(t)*bw*W(t) - mu*I(t),
diff(S(t), t) = -bi*S(t)*I(t) - S(t)*mu - S(t)*bw*W(t) + R(t)*a + mu,
diff(W(t), t) = xi*I(t) - xi*W(t),
y2(t) = S(t) + R(t) + I(t),
y(t) = k*I(t)
];
IdentifiabilityODE(sigma, GetParameters(sigma));