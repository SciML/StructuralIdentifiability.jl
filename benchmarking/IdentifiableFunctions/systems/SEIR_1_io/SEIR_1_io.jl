# SEIR_1_io
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	R'(t) = -I(t)*gamma*psi + I(t)*gamma + gamma*Q(t),
	Q'(t) = I(t)*psi - gamma*Q(t),
	I'(t) = I(t)*gamma*psi - I(t)*gamma - I(t)*psi + v*E(t),
	E'(t) = beta*I(t)*S(t) - v*E(t),
	S'(t) = -beta*I(t)*S(t),
	y1(t) = Q(t)
)
