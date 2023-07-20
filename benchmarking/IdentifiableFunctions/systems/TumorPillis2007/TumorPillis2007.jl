# TumorPillis2007
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	C'(t) = alpha2 - KC*M(t)*C(t) - beta*C(t),
	N'(t) = (-p*h*N(t)*T(t) - p*N(t)*T(t)^2 - h*KN*N(t)*M(t) - h*N(t)*f + h*alpha1 - KN*N(t)*M(t)*T(t) + N(t)*g*T(t) - N(t)*T(t)*f + alpha1*T(t))//(h + T(t)),
	I'(t) = (u2(t)*gt + u2(t)*T(t) + w*gt*I(t)*L(t) + w*T(t)*I(t)*L(t) - muI*gt*I(t) - muI*T(t)*I(t) + pt*T(t)*L(t))//(gt + T(t)),
	M'(t) = u1(t) - gamma*M(t),
	T'(t) = -b*T(t)^2*a - D(t)*T(t) - KT*M(t)*T(t) - c1*N(t)*T(t) + T(t)*a,
	L'(t) = (u1(t)*gI + u1(t)*I(t) + r2*gI*T(t)*C(t) + r2*T(t)*I(t)*C(t) + pI*I(t)*L(t) - KL*M(t)*gI*L(t) - KL*M(t)*I(t)*L(t) - gI*q*T(t)*L(t) - gI*m*L(t) - gI*L(t)^2*ucte - q*T(t)*I(t)*L(t) - m*I(t)*L(t) - I(t)*L(t)^2*ucte)//(gI + I(t)),
	y2(t) = N(t),
	y1(t) = L(t),
	y3(t) = M(t)
)
