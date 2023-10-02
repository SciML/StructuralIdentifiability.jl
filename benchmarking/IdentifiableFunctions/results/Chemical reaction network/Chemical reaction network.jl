# Chemical reaction network
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	x1'(t) = k4*x6(t) + k2*x4(t) - k1*x1(t)*x2(t),
	x2'(t) = k3*x4(t) + k2*x4(t) - k1*x1(t)*x2(t),
	x3'(t) = k5*x6(t) + k3*x4(t) - k6*x5(t)*x3(t),
	x4'(t) = -k3*x4(t) - k2*x4(t) + k1*x1(t)*x2(t),
	x5'(t) = k5*x6(t) + k4*x6(t) - k6*x5(t)*x3(t),
	x6'(t) = -k5*x6(t) - k4*x6(t) + k6*x5(t)*x3(t),
	y1(t) = x3(t),
	y2(t) = x2(t)
)
