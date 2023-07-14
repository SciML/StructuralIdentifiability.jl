# QWWC
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	z'(t) = -c*z(t) - w(t)*d + x(t)*y(t),
	x'(t) = -x(t)*a + z(t)*y(t) + a*y(t),
	w'(t) = e*z(t) - w(t)*f + x(t)*y(t),
	y'(t) = b*x(t) + b*y(t) - x(t)*z(t),
	g(t) = x(t)
)
