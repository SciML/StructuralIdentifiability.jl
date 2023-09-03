# Crauste_SI
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	M'(t) = S(t)*delta_LM - M(t)*mu_M,
	P'(t) = rho_P*P(t)^2 - S(t)*P(t)*mu_PL - E(t)*mu_PE*P(t) - P(t)*mu_P,
	E'(t) = -mu_EE*E(t)^2 + delta_NE*N(t)*P(t) - delta_EL*E(t) + E(t)*P(t)*rho_E,
	N'(t) = -delta_NE*N(t)*P(t) - mu_N*N(t),
	S'(t) = -mu_LE*S(t)*E(t) + delta_EL*S(t) - S(t)^2*mu_LL - S(t)*delta_LM,
	y3(t) = M(t),
	y1(t) = N(t),
	y2(t) = S(t) + E(t)
)
