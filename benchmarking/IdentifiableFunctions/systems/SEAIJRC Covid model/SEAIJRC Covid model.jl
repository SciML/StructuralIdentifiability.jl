# SEAIJRC Covid model
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	A'(t) = -E(t)*k*r + E(t)*k - A(t)*g1,
	Ninv'(t) = 0,
	I'(t) = -alpha*I(t) + E(t)*k*r - g1*I(t),
	C'(t) = alpha*I(t),
	J'(t) = alpha*I(t) - g2*J(t),
	S'(t) = -b*S(t)*Ninv(t)*A(t)*q - b*S(t)*Ninv(t)*I(t) - b*S(t)*Ninv(t)*J(t),
	E'(t) = b*S(t)*Ninv(t)*A(t)*q + b*S(t)*Ninv(t)*I(t) + b*S(t)*Ninv(t)*J(t) - E(t)*k,
	y(t) = C(t),
	y2(t) = Ninv(t)
)
