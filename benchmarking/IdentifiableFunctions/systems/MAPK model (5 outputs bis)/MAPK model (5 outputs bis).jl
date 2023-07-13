# MAPK model (5 outputs bis)
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	KS00'(t) = FS11(t)*gamma1100 + gamma1000*FS10(t) + FS01(t)*gamma0100 + b00*KS00(t) - S00(t)*K(t)*a00,
	S01'(t) = -FS11(t)*beta11 - FS11(t)*gamma1100 - FS11(t)*gamma1101 - FS11(t)*gamma1110 + F(t)*alpha11*S11(t),
	KS01'(t) = c0001*KS00(t) + FS11(t)*gamma1101 - F(t)*alpha01*S01(t) + FS01(t)*beta01 - K(t)*S01(t)*a01 + KS01(t)*b01,
	FS11'(t) = -c0111*KS01(t) + K(t)*S01(t)*a01 - KS01(t)*b01,
	KS10'(t) = -a10*K(t)*S10(t) + FS11(t)*gamma1110 - F(t)*alpha10*S10(t) + KS10(t)*b10 + beta10*FS10(t) + c0010*KS00(t),
	FS01'(t) = FS11(t)*beta11 - F(t)*alpha11*S11(t) + KS10(t)*c1011 + c0111*KS01(t) + c0011*KS00(t),
	S00'(t) = F(t)*alpha10*S10(t) - gamma1000*FS10(t) - beta10*FS10(t),
	K'(t) = a10*K(t)*S10(t) - KS10(t)*b10 - KS10(t)*c1011,
	FS10'(t) = -c0001*KS00(t) - b00*KS00(t) + S00(t)*K(t)*a00 - c0011*KS00(t) - c0010*KS00(t),
	S10'(t) = c0001*KS00(t) - a10*K(t)*S10(t) + KS10(t)*b10 + KS10(t)*c1011 + b00*KS00(t) + c0111*KS01(t) - S00(t)*K(t)*a00 - K(t)*S01(t)*a01 + c0011*KS00(t) + c0010*KS00(t) + KS01(t)*b01,
	S11'(t) = FS11(t)*beta11 + FS11(t)*gamma1100 + FS11(t)*gamma1101 + FS11(t)*gamma1110 - F(t)*alpha10*S10(t) - F(t)*alpha11*S11(t) - F(t)*alpha01*S01(t) + gamma1000*FS10(t) + FS01(t)*beta01 + FS01(t)*gamma0100 + beta10*FS10(t),
	F'(t) = F(t)*alpha01*S01(t) - FS01(t)*beta01 - FS01(t)*gamma0100,
	y4(t) = S11(t),
	y2(t) = S00(t),
	y1(t) = F(t),
	y3(t) = S01(t) + S10(t),
	y0(t) = K(t)
)
