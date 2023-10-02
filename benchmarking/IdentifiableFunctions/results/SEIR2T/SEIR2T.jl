# SEIR2T
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	S'(t) = (-b*In(t)*S(t))//N(t),
	E'(t) = (b*In(t)*S(t) - N(t)*nu*E(t))//N(t),
	In'(t) = -a*In(t) + nu*E(t),
	N'(t) = 0,
	Cu'(t) = nu*E(t),
	y1(t) = Cu(t),
	y2(t) = N(t)
)
