# SIR 21
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	I'(t) = (S(t)*beta*I(t) - N(t)*mu*I(t) - N(t)*r*I(t))//N(t),
	N'(t) = 0,
	D'(t) = mu*I(t),
	C'(t) = S(t)*pp - q*C(t),
	S'(t) = (-S(t)*N(t)*pp - S(t)*beta*I(t) + N(t)*q*C(t))//N(t),
	R'(t) = r*I(t),
	y1(t) = N(t),
	y3(t) = C(t),
	y2(t) = D(t)
)
