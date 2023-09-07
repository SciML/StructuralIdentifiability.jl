# SEIRT
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	S'(t) = (-beta*I(t)*S(t))//N(t),
	E'(t) = (-alpha*N(t)*E(t) + beta*I(t)*S(t))//N(t),
	I'(t) = -lambda*I(t) + alpha*E(t),
	R'(t) = lambda*I(t),
	N'(t) = 0,
	y1(t) = I(t),
	y2(t) = N(t)
)
