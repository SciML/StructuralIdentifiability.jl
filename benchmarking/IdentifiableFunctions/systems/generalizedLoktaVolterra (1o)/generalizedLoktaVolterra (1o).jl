# generalizedLoktaVolterra (1o)
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	x1'(t) = r1*x1(t) + beta11*x1(t)^2 + x2(t)*beta12*x1(t),
	x2'(t) = beta22*x2(t)^2 + beta21*x2(t)*x1(t) + x2(t)*r2,
	y1(t) = x1(t)
)
