# CD8 T cell differentiation
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	N'(t) = -delta_NE*N(t)*P(t) - mu_N*N(t),
	E'(t) = -mu_EE*E(t)^2 + delta_NE*N(t)*P(t) - delta_EL*E(t) + E(t)*P(t)*rho_E,
	S'(t) = -mu_LE*S(t)*E(t) + delta_EL*S(t) - S(t)^2*mu_LL - S(t)*delta_LM,
	M'(t) = S(t)*delta_LM - M(t)*mu_M,
	P'(t) = rho_P*P(t)^2 - S(t)*P(t)*mu_PL - E(t)*mu_PE*P(t) - P(t)*mu_P,
	y1(t) = N(t),
	y2(t) = S(t) + E(t),
	y3(t) = M(t)
)
