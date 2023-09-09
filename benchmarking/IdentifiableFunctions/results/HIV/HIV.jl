# HIV
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	x'(t) = lm - x(t)*d - x(t)*v(t)*beta,
	y'(t) = x(t)*v(t)*beta - a*y(t),
	v'(t) = k*y(t) - v(t)*u,
	w'(t) = -b*w(t) + c*w(t)*x(t)*y(t) - c*w(t)*q*y(t),
	z'(t) = c*w(t)*q*y(t) - h*z(t),
	y1(t) = w(t),
	y2(t) = z(t)
)
