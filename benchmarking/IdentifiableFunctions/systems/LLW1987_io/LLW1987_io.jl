# LLW1987_io
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	x2'(t) = -p3*x2(t) + p4*u(t),
	x1'(t) = p2*u(t) - p1*x1(t),
	x3'(t) = p2*x2(t)*u(t) - p3*x3(t) + p4*u(t)*x1(t) - x3(t)*p1,
	y1(t) = x3(t)
)
