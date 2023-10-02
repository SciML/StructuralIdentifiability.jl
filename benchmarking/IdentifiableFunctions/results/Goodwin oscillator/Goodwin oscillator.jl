# Goodwin oscillator
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	x1'(t) = (-b*c*x1(t) - b*x1(t)*x4(t) + 1)//(c + x4(t)),
	x2'(t) = alpha*x1(t) - beta*x2(t),
	x3'(t) = gama*x2(t) - delta*x3(t),
	x4'(t) = (gama*sigma*x2(t)*x4(t) - delta*sigma*x3(t)*x4(t))//x3(t),
	y(t) = x1(t)
)
