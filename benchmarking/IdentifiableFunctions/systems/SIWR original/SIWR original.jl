# SIWR original
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	I'(t) = bi*S(t)*I(t) - gam*I(t) + S(t)*bw*W(t) - mu*I(t),
	R'(t) = gam*I(t) - R(t)*mu - R(t)*a,
	S'(t) = -bi*S(t)*I(t) - S(t)*mu - S(t)*bw*W(t) + R(t)*a + mu,
	W'(t) = xi*I(t) - xi*W(t),
	y(t) = k*I(t)
)
