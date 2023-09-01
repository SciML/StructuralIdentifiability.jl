# SEUIR
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	R'(t) = d*U(t) + d*I(t),
	U'(t) = -w*E(t) + E(t)*z - d*U(t),
	E'(t) = (S(t)*U(t)*beta + S(t)*beta*I(t) - E(t)*N*z)//N,
	S'(t) = (-S(t)*U(t)*beta - S(t)*beta*I(t))//N,
	I'(t) = w*E(t) - d*I(t),
	y1(t) = I(t)
)
