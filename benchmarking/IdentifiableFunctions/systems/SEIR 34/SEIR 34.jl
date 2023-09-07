# SEIR 34
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	A'(t) = 0,
	N'(t) = 0,
	S'(t) = (-S(t)*N(t)*mu - S(t)*r*beta*I(t) + N(t)*A(t))//N(t),
	E'(t) = (S(t)*r*beta*I(t) - epsilon*E(t)*N(t) - E(t)*N(t)*mu)//N(t),
	I'(t) = -gamma*I(t) + epsilon*E(t) - mu*I(t),
	R'(t) = gamma*I(t) - R(t)*mu,
	y1(t) = K*I(t),
	y2(t) = A(t),
	y3(t) = N(t)
)
