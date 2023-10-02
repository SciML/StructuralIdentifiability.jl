# SIR 24
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	A'(t) = 0,
	S'(t) = (-c*S(t)*phi*I(t) + S(t)^2 + S(t)*A(t) - S(t)*mu + S(t)*I(t) + A(t)*I(t) - mu*I(t))//(S(t) + I(t)),
	I'(t) = (-u1(t)*I(t) - gamma*S(t)*I(t) - gamma*I(t)^2 + c*S(t)*phi*I(t) - S(t)*mu*I(t) - mu*I(t)^2)//(S(t) + I(t)),
	R'(t) = (u1(t)*I(t) + gamma*S(t)*I(t) + gamma*I(t)^2 - S(t)*R(t)*mu - R(t)*mu*I(t))//(S(t) + I(t)),
	y1(t) = K*I(t)
)
