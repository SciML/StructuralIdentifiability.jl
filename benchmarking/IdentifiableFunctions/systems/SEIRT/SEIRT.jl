# SEIRT
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	I'(t) = -lambda*I(t) + alpha*E(t),
	E'(t) = (-alpha*N(t)*E(t) + beta*I(t)*S(t))//N(t),
	N'(t) = 0,
	R'(t) = lambda*I(t),
	S'(t) = (-beta*I(t)*S(t))//N(t),
	y2(t) = N(t),
	y1(t) = I(t)
)
