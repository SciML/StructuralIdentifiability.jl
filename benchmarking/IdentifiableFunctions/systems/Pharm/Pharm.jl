# Pharm
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	x0'(t) = (-ka*n*x0(t) - ka*kc*a1*x0(t) + ka*kc*a1*x1(t) - ka*a1*x0(t)^2 + ka*a1*x0(t)*x1(t) - kc*a1*x0(t)*x2(t) + kc*a1*x1(t)*x2(t))//(ka*kc + ka*x0(t) + kc*x2(t)),
	x1'(t) = a2*x0(t) - a2*x1(t),
	x3'(t) = -b2*x3(t) + b2*x2(t),
	x2'(t) = (ka*kc*b1*x3(t) - ka*kc*b1*x2(t) + ka*b1*x0(t)*x3(t) - ka*b1*x0(t)*x2(t) - n*kc*x2(t) + kc*b1*x3(t)*x2(t) - kc*b1*x2(t)^2)//(ka*kc + ka*x0(t) + kc*x2(t)),
	y1(t) = x0(t)
)
