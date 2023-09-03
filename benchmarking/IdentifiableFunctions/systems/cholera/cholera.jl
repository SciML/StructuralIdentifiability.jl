# cholera
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	R'(t) = gam*I(t) - R(t)*mu - R(t)*a,
	I'(t) = bi*S(t)*I(t) - gam*I(t) + S(t)*bw*W(t) - mu*I(t),
	S'(t) = -bi*S(t)*I(t) - S(t)*mu - S(t)*bw*W(t) + R(t)*a + mu,
	W'(t) = xi*I(t) - xi*W(t),
	y2(t) = S(t) + R(t) + I(t),
	y1(t) = k*I(t)
)
