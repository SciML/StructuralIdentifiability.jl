# Bilirubin2_io
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	x1'(t) = -k01*x1(t) - k31*x1(t) - k21*x1(t) + k12*x2(t) + x3(t)*k13 + u(t) - x1(t)*k41 + x4(t)*k14,
	x2'(t) = k21*x1(t) - k12*x2(t),
	x3'(t) = k31*x1(t) - x3(t)*k13,
	x4'(t) = x1(t)*k41 - x4(t)*k14,
	y1(t) = x1(t)
)
