# Pivastatin
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	x1'(t) = -r3*x1(t) + k3*x3(t) + r1*x2(t) - T0*k1*x1(t) + k1*x1(t)*x2(t),
	x2'(t) = -r1*x2(t) - k2*x2(t) + T0*k1*x1(t) - k1*x1(t)*x2(t),
	x3'(t) = r3*x1(t) - k3*x3(t) - k4*x3(t) + k2*x2(t),
	y1(t) = k*x3(t) + k*x2(t)
)
