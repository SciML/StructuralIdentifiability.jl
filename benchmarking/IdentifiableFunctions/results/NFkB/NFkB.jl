# NFkB
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	x8'(t) = 1//2*x9(t) - c5*x8(t),
	x3'(t) = x2(t)*k2*u(t)*x8(t) + x2(t)*k3 - x3(t)*k_deg,
	x10'(t) = -1//5*x10(t)*x2(t) - x10(t)*i1a - 1//2*x10(t)*x6(t) - 1//10000*x10(t) + 1//2000*x11(t) + x12(t)*c4a,
	x1'(t) = -x1(t)*k_deg - x1(t)*u(t)*k1 + k_prod,
	x12'(t) = 1//2000000*x7(t) - c3a*x12(t),
	x9'(t) = -1//2500*x9(t) - 1//2000000*x7(t),
	x6'(t) = -1//2*x10(t)*x6(t) + 1//50000*x13(t) - i1*x6(t) + x5(t)*t2,
	x4'(t) = 1//5*x10(t)*x2(t) - x4(t)*t1,
	x7'(t) = 5*i1*x6(t) - 1//2*x11(t)*x7(t),
	x11'(t) = 5*x10(t)*i1a - 1//2*x11(t)*x7(t) - 1//400*x11(t),
	x5'(t) = x2(t)*x13(t) - x5(t)*t2,
	x14'(t) = -5*e2a*x14(t) + 1//2*x11(t)*x7(t),
	x15'(t) = -1//2500*x15(t) + 1//2000000*x7(t),
	x13'(t) = e2a*x14(t) + 1//2*x10(t)*x6(t) - x2(t)*x13(t) - 1//50000*x13(t),
	x2'(t) = -1//5*x10(t)*x2(t) - x2(t)*x13(t) - x2(t)*k_deg - x2(t)*k2*u(t)*x8(t) - x2(t)*k3 + x1(t)*u(t)*k1 + x4(t)*t1 + x5(t)*t2,
	y3(t) = x9(t),
	y6(t) = x12(t),
	y2(t) = x10(t) + x13(t),
	y4(t) = x2(t) + x3(t) + x1(t),
	y5(t) = x2(t),
	y1(t) = x7(t)
)
