# SIR 6
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	N'(t) = 0,
	S'(t) = (-beta*I(t)*S(t))//N(t),
	I'(t) = (-N(t)*I(t)*gamma + beta*I(t)*S(t))//N(t),
	R'(t) = I(t)*gamma,
	y1(t) = I(t)*K,
	y2(t) = N(t)
)
