# SLIQR
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	S'(t) = -b*In(t)*S(t)*Ninv - S(t)*Ninv*u(t),
	In'(t) = -In(t)*g + s*Q(t) + a*L(t),
	L'(t) = b*In(t)*S(t)*Ninv - a*L(t),
	Q'(t) = -e*In(t)*g + In(t)*g - s*Q(t),
	y(t) = In(t)*Ninv
)
