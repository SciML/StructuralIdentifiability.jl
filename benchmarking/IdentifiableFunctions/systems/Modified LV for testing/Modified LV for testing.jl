# Modified LV for testing
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	x1'(t) = a*x1(t) + b*x1(t) - x2(t)*c*x1(t),
	x2'(t) = -a*b*x2(t) + d*x2(t)*x1(t),
	y1(t) = x1(t)
)
