# St
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	S'(t) = -e*S(t) - S(t)*d*W(t) + S(t)*r - S(t)*a*W(t) + R(t)*g,
	R'(t) = e*S(t) + rR*R(t) + S(t)*a*W(t) - dr*R(t)*W(t) - R(t)*g,
	W'(t) = T*Dd - W(t)*Dd,
	y1(t) = S(t) + R(t),
	y2(t) = T
)
