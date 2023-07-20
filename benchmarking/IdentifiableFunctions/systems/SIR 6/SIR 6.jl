# SIR 6
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	I'(t) = (-N(t)*I(t)*gamma + beta*I(t)*S(t))//N(t),
	R'(t) = I(t)*gamma,
	S'(t) = (-beta*I(t)*S(t))//N(t),
	N'(t) = 0,
	y2(t) = N(t),
	y1(t) = I(t)*K
)
