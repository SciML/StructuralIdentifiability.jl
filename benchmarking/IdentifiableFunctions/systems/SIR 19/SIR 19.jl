# SIR 19
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	N'(t) = 0,
	S'(t) = (-S(t)*N(t)*pp - S(t)*beta*I(t) + N(t)*q*C(t))//N(t),
	I'(t) = (S(t)*beta*I(t) - N(t)*mu*I(t) - N(t)*r*I(t))//N(t),
	R'(t) = r*I(t),
	C'(t) = S(t)*pp - q*C(t),
	D'(t) = mu*I(t),
	y1(t) = N(t),
	y2(t) = D(t),
	y3(t) = C(t)
)
