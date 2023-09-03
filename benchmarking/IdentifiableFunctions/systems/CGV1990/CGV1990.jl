# CGV1990
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	q35'(t) = k5*R*q3(t)*V3 - k5*q3(t)*q35(t) - k6*q35(t),
	q3'(t) = (5*k5*q36(t)*V36*q3(t) - 5*k5*S*V36^2*q3(t) - k5*R*q3(t)*V3^2 + k5*q3(t)*q35(t)*V3 + k3*q1(t)*V3 + q36(t)*k6*V3 - k4*q3(t)*V3 + k6*q35(t)*V3)//V3,
	q1'(t) = -k3*q1(t) - q1(t)*k7 + k4*q3(t) + u(t),
	q7'(t) = q1(t)*k7,
	q36'(t) = (-5*k5*q36(t)*V36*q3(t) + 5*k5*S*V36^2*q3(t) - q36(t)*k6*V3)//V3,
	y1(t) = q7(t)
)
