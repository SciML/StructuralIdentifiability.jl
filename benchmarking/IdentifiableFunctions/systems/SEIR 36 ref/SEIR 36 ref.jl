# SEIR 36 ref
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	N'(t) = 0,
	nu'(t) = 0,
	q'(t) = 0,
	S'(t) = (nu(t)*N(t)^2 - S(t)*mu_0*N(t) - S(t)*q(t)*Di(t)*beta_d - S(t)*beta*I(t))//N(t),
	E'(t) = (S(t)*q(t)*Di(t)*beta_d + S(t)*beta*I(t) - mu_0*E(t)*N(t) - E(t)*N(t)*s - E(t)*N(t)*phi_e)//N(t),
	I'(t) = -gamma*I(t) - phi*I(t) - mu_0*I(t) + E(t)*s - mu_i*I(t),
	De'(t) = -De(t)*mu_0 - De(t)*s_d + E(t)*phi_e,
	Di'(t) = De(t)*s_d + phi*I(t) - mu_0*Di(t) - mu_d*Di(t) - Di(t)*gamma_d,
	R'(t) = gamma*I(t) - mu_0*R(t) + Di(t)*gamma_d,
	F'(t) = mu_d*Di(t) + mu_i*I(t),
	y1(t) = De(t),
	y2(t) = Di(t),
	y5(t) = F(t),
	y3(t) = N(t),
	y4(t) = nu(t),
	y6(t) = q(t)
)
