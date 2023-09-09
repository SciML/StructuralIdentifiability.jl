# SIWR with extra output
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	S'(t) = -bi*S(t)*I(t) - S(t)*mu - S(t)*bw*W(t) + R(t)*a + mu,
	I'(t) = bi*S(t)*I(t) - gam*I(t) + S(t)*bw*W(t) - mu*I(t),
	W'(t) = xi*I(t) - xi*W(t),
	R'(t) = gam*I(t) - R(t)*mu - R(t)*a,
	y(t) = k*I(t),
	y2(t) = S(t) + R(t) + I(t)
)
