# SIRS forced
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	s'(t) = -b1*b0*s(t)*x1(t)*i(t) - b0*s(t)*i(t) - s(t)*mu + mu + g*r(t),
	i'(t) = -nu*i(t) + b1*b0*s(t)*x1(t)*i(t) + b0*s(t)*i(t) - mu*i(t),
	r'(t) = nu*i(t) - mu*r(t) - g*r(t),
	x1'(t) = -M*x2(t),
	x2'(t) = M*x1(t),
	y1(t) = i(t),
	y2(t) = r(t)
)
