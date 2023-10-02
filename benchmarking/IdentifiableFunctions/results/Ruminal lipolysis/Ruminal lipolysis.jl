# Ruminal lipolysis
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	x1'(t) = (-x5(t)*x1(t))//(k2 + x1(t)),
	x2'(t) = (-k4*k2*x2(t) - k4*x1(t)*x2(t) + 2//3*x5(t)*x1(t))//(k2 + x1(t)),
	x3'(t) = -k4*x3(t) + 1//2*k4*x2(t),
	x4'(t) = (k4*k2*x3(t) + 1//2*k4*k2*x2(t) + k4*x3(t)*x1(t) + 1//2*k4*x1(t)*x2(t) + 1//3*x5(t)*x1(t))//(k2 + x1(t)),
	x5'(t) = -k3*x5(t),
	y1(t) = x1(t),
	y2(t) = x3(t) + x2(t),
	y3(t) = x4(t)
)
