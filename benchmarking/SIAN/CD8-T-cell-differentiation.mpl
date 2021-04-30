S(t) + E(t)diff(N(t), t) = -delta_NE*N(t)*P(t) - mu_N*N(t),
diff(S(t), t) = -mu_LE*S(t)*E(t) + delta_EL*S(t) - S(t)^2*mu_LL - S(t)*delta_LM,
diff(M(t), t) = S(t)*delta_LM - M(t)*mu_M,
diff(P(t), t) = rho_P*P(t)^2 - S(t)*P(t)*mu_PL - E(t)*mu_PE*P(t) - P(t)*mu_P,
diff(E(t), t) = -mu_EE*E(t)^2 + delta_NE*N(t)*P(t) - delta_EL*E(t) + E(t)*P(t)*rho_E,
y1(t) = N(t),
y3(t) = M(t),
y2(t) = S(t) + E(t)
];
IdentifiabilityODE(sigma, GetParameters(sigma));