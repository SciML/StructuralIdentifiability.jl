function model = SMTH()
syms KS01 KS10 FS01 FS11 K FS10 KS00 S10 S00 S01 S11 F
syms c0001 a10 gamma1000 alpha10 b00 beta11 c0111 beta10 alpha11 beta01 alpha01 gamma1100 c0011 c0010 b10 gamma1101 a00 b01 c1011 a01 gamma1110 gamma0100
syms KS010 KS100 FS010 FS110 K0 FS100 KS000 S100 S000 S010 S110 F0
model.sym.p = [c0001; a10; gamma1000; alpha10; b00; beta11; c0111; beta10; alpha11; beta01; alpha01; gamma1100; c0011; c0010; b10; gamma1101; a00; b01; c1011; a01; gamma1110; gamma0100; KS010; KS100; FS010; FS110; K0; FS100; KS000; S100; S000; S010; S110; F0];
model.sym.x = [KS01; KS10; FS01; FS11; K; FS10; KS00; S10; S00; S01; S11; F];
model.sym.g = [];
model.sym.x0 = [KS010; KS100; FS010; FS110; K0; FS100; KS000; S100; S000; S010; S110; F0];
model.sym.xdot = [c0001*KS00 + FS11*gamma1101 - F*alpha01*S01 + FS01*beta01 - K*S01*a01 + KS01*b01
-a10*K*S10 + FS11*gamma1110 - F*alpha10*S10 + KS10*b10 + beta10*FS10 + c0010*KS00
FS11*beta11 - F*alpha11*S11 + KS10*c1011 + c0111*KS01 + c0011*KS00
-c0111*KS01 + K*S01*a01 - KS01*b01
a10*K*S10 - KS10*b10 - KS10*c1011
-c0001*KS00 - b00*KS00 + S00*K*a00 - c0011*KS00 - c0010*KS00
FS11*gamma1100 + gamma1000*FS10 + FS01*gamma0100 + b00*KS00 - S00*K*a00
c0001*KS00 - a10*K*S10 + KS10*b10 + KS10*c1011 + b00*KS00 + c0111*KS01 - S00*K*a00 - K*S01*a01 + c0011*KS00 + c0010*KS00 + KS01*b01
F*alpha10*S10 - gamma1000*FS10 - beta10*FS10
-FS11*beta11 - FS11*gamma1100 - FS11*gamma1101 - FS11*gamma1110 + F*alpha11*S11
FS11*beta11 + FS11*gamma1100 + FS11*gamma1101 + FS11*gamma1110 - F*alpha10*S10 - F*alpha11*S11 - F*alpha01*S01 + gamma1000*FS10 + FS01*beta01 + FS01*gamma0100 + beta10*FS10
F*alpha01*S01 - FS01*beta01 - FS01*gamma0100];
model.sym.y = [F
S10
S11
S00
S01];
end