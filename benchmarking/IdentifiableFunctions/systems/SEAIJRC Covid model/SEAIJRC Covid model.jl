# SEAIJRC Covid model
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	S'(t) = -b*S(t)*Ninv(t)*A(t)*q - b*S(t)*Ninv(t)*I(t) - b*S(t)*Ninv(t)*J(t),
	E'(t) = b*S(t)*Ninv(t)*A(t)*q + b*S(t)*Ninv(t)*I(t) + b*S(t)*Ninv(t)*J(t) - E(t)*k,
	A'(t) = -E(t)*k*r + E(t)*k - A(t)*g1,
	I'(t) = -alpha*I(t) + E(t)*k*r - g1*I(t),
	J'(t) = alpha*I(t) - g2*J(t),
	C'(t) = alpha*I(t),
	Ninv'(t) = 0,
	y(t) = C(t),
	y2(t) = Ninv(t)
)
