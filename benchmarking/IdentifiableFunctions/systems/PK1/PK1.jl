# PK1
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	x2'(t) = k5*x4(t) - k3*x2(t) - k6*x2(t) + k1*x1(t) - k7*x2(t),
	x1'(t) = u1(t) - k2*x1(t) - k1*x1(t),
	x3'(t) = k3*x2(t) - k4*x3(t) + k2*x1(t),
	x4'(t) = -k5*x4(t) + k6*x2(t),
	y2(t) = s3*x3(t),
	y1(t) = s2*x2(t)
)
