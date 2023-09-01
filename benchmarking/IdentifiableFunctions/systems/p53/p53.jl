# p53
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	x3'(t) = (u1(t)*p18*p16*x3(t)*x1(t) - p16*x3(t)*x1(t) - p17*x3(t)*p15 + p17*p14 - x3(t)^2*p15 + x3(t)*p14)//(p17 + x3(t)),
	x1'(t) = (-u1(t)*p4*p6*x1(t)^2 - u1(t)*p4*x1(t)^2 + u1(t)*p1*x1(t)*x4(t) + u1(t)*p1*p5*x4(t) - u1(t)*p3*x1(t)^2 - u1(t)*p3*x1(t)*p5 - p7*p4*x1(t)^2 + p7*p1*x1(t)*x4(t) + p7*p1*p5*x4(t) - p7*p3*x1(t)^2 - p7*p3*x1(t)*p5)//(u1(t)*x1(t) + u1(t)*p5 + p7*x1(t) + p7*p5),
	x4'(t) = (u1(t)*p23*p25*p22^4*x1(t)*x2(t)*p24 - u1(t)*p23*p25*p22^4*x2(t) + u1(t)*p23*p25*x1(t)*x2(t)*p24 - u1(t)*p23*p25*x2(t) - u1(t)*p23*p22^4*x1(t)*p24 + u1(t)*p23*p22^4 - u1(t)*p23*x1(t)*p24 + u1(t)*p23 - p25*p22^4*p21*p24 + p25*p22^4*p21 + p25*p22^4*x1(t)*x2(t)*p24 - p25*p22^4*x2(t) - p25*p21*x3(t)^4*p24 + p25*p21*x3(t)^4 + p25*x1(t)*x2(t)*p24 - p25*x2(t) - p20*p22^8*x4(t) + p20*p22^8 - p20*p22^4*x3(t)^4*x4(t) + p20*p22^4*x3(t)^4 - p20*p22^4*x4(t) + p20*p22^4 - p20*x3(t)^4*x4(t) + p20*x3(t)^4 + p22^8*p21*x3(t)^4 + p22^4*p21*x3(t)^8 + p22^4*p21*x3(t)^4 + p22^4*p21*p24 - p22^4*p21 - p22^4*x1(t)*p24 + p22^4 + p21*x3(t)^8 + p21*x3(t)^4*p24 - p21*x3(t)^4 - x1(t)*p24 + 1)//(p22^8 + p22^4*x3(t)^4 + p22^4 + x3(t)^4),
	x2'(t) = (u1(t)*p8*p11 + u1(t)*p8*x2(t) - u1(t)*p10*p12*x1(t)*x2(t) - u1(t)*p10*x1(t)*x2(t) - u1(t)*p11*p9*x2(t) - u1(t)*p9*x2(t)^2 + p13*p8*p11 + p13*p8*x2(t) - p13*p10*x1(t)*x2(t) - p13*p11*p9*x2(t) - p13*p9*x2(t)^2)//(u1(t)*p11 + u1(t)*x2(t) + p13*p11 + p13*x2(t)),
	y2(t) = x2(t),
	y3(t) = x3(t),
	y4(t) = x4(t),
	y1(t) = x1(t)
)
