# HIV2_io
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	x1'(t) = -b*x1(t)*x4(t) - d*x1(t) + s,
	x2'(t) = b*q1*x1(t)*x4(t) - w1*x2(t) - k1*x2(t),
	x3'(t) = b*x1(t)*q2*x4(t) - x3(t)*w2 + k1*x2(t),
	x4'(t) = -c*x4(t) + k2*x3(t),
	y1(t) = x1(t),
	y2(t) = x4(t)
)
