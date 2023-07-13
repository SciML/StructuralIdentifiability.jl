# HIV
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	w'(t) = -b*w(t) + c*w(t)*x(t)*y(t) - c*w(t)*q*y(t),
	v'(t) = k*y(t) - v(t)*u,
	x'(t) = lm - x(t)*d - x(t)*v(t)*beta,
	z'(t) = c*w(t)*q*y(t) - h*z(t),
	y'(t) = x(t)*v(t)*beta - a*y(t),
	y2(t) = z(t),
	y1(t) = w(t)
)
