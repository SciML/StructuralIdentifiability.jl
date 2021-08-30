with(DifferentialAlgebra):
ring_diff := DifferentialRing(blocks = [[KS00, KS10, S01, FS01, FS11, K, KS01, FS10, S11, S00, S10, F], [y0, y4, y2, y1, y3] ], derivations = [t]):
sys := [
diff(KS00(t), t) - (FS11(t)*gamma1100 + gamma1000*FS10(t) + FS01(t)*gamma0100 + b00*KS00(t) - S00(t)*K(t)*a00),
diff(KS10(t), t) - (-a10*K(t)*S10(t) + FS11(t)*gamma1110 - F(t)*alpha10*S10(t) + KS10(t)*b10 + beta10*FS10(t) + c0010*KS00(t)),
diff(S01(t), t) - (-FS11(t)*beta11 - FS11(t)*gamma1100 - FS11(t)*gamma1101 - FS11(t)*gamma1110 + F(t)*alpha11*S11(t)),
diff(FS01(t), t) - (FS11(t)*beta11 - F(t)*alpha11*S11(t) + KS10(t)*c1011 + c0111*KS01(t) + c0011*KS00(t)),
diff(FS11(t), t) - (-c0111*KS01(t) + K(t)*S01(t)*a01 - KS01(t)*b01),
diff(K(t), t) - (a10*K(t)*S10(t) - KS10(t)*b10 - KS10(t)*c1011),
diff(KS01(t), t) - (c0001*KS00(t) + FS11(t)*gamma1101 - F(t)*alpha01*S01(t) + FS01(t)*beta01 - K(t)*S01(t)*a01 + KS01(t)*b01),
diff(FS10(t), t) - (-c0001*KS00(t) - b00*KS00(t) + S00(t)*K(t)*a00 - c0011*KS00(t) - c0010*KS00(t)),
diff(S11(t), t) - (FS11(t)*beta11 + FS11(t)*gamma1100 + FS11(t)*gamma1101 + FS11(t)*gamma1110 - F(t)*alpha10*S10(t) - F(t)*alpha11*S11(t) - F(t)*alpha01*S01(t) + gamma1000*FS10(t) + FS01(t)*beta01 + FS01(t)*gamma0100 + beta10*FS10(t)),
diff(S00(t), t) - (F(t)*alpha10*S10(t) - gamma1000*FS10(t) - beta10*FS10(t)),
diff(S10(t), t) - (c0001*KS00(t) - a10*K(t)*S10(t) + KS10(t)*b10 + KS10(t)*c1011 + b00*KS00(t) + c0111*KS01(t) - S00(t)*K(t)*a00 - K(t)*S01(t)*a01 + c0011*KS00(t) + c0010*KS00(t) + KS01(t)*b01),
diff(F(t), t) - (F(t)*alpha01*S01(t) - FS01(t)*beta01 - FS01(t)*gamma0100),
y0(t) - (K(t)),
y4(t) - (S11(t)),
y2(t) - (S00(t)),
y1(t) - (F(t)),
y3(t) - (S01(t) + S10(t))
];
res := CodeTools[CPUTime](RosenfeldGroebner(sys, ring_diff, singsol=none));