# SEUIR
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	S'(t) = (-S(t)*U(t)*beta - S(t)*beta*I(t))//N,
	E'(t) = (S(t)*U(t)*beta + S(t)*beta*I(t) - E(t)*N*z)//N,
	U'(t) = -w*E(t) + E(t)*z - d*U(t),
	I'(t) = w*E(t) - d*I(t),
	R'(t) = d*U(t) + d*I(t),
	y1(t) = I(t)
)
