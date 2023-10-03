# Biohydrogenation_io
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	x4'(t) = (-k5*x4(t))//(k6 + x4(t)),
	x5'(t) = (k5*k8*x4(t) + k5*x6(t)*x4(t) + k5*x5(t)*x4(t) - k6*x5(t)*k7 - x5(t)*k7*x4(t))//(k8*k6 + k8*x4(t) + k6*x6(t) + k6*x5(t) + x6(t)*x4(t) + x5(t)*x4(t)),
	x6'(t) = (-k8*k9*k10*x6(t) + k8*k9*x6(t)^2 - k9*k10*x6(t)^2 - k9*k10*x6(t)*x5(t) + k9*x6(t)^3 + k9*x6(t)^2*x5(t) + k10*x5(t)*k7)//(k8*k10 + k10*x6(t) + k10*x5(t)),
	x7'(t) = (k9*k10*x6(t) - k9*x6(t)^2)//k10,
	y1(t) = x4(t),
	y2(t) = x5(t)
)
