# JAK-STAT 1
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	x1'(t) = t6*x2(t) - t5*x1(t) - 2*u(t)*x1(t)*t1,
	x2'(t) = -t6*x2(t) + t5*x1(t),
	x3'(t) = x6(t)*x3(t)*t2 - 3*x3(t)*t2 + 2*u(t)*x1(t)*t1,
	x4'(t) = -t3*x4(t) - x6(t)*x3(t)*t2 + 3*x3(t)*t2,
	x5'(t) = t3*x4(t) - x5(t)*t4,
	x6'(t) = (-x6(t)*x3(t)*x10(t)*t7*t13 - x6(t)*x3(t)*t7 - 92*x6(t)*x10(t)*x1(t)*t8*t13^2 - 92*x6(t)*x10(t)*t8*t13 - 92*x6(t)*x1(t)*t8*t13 - x6(t)*x1(t)*t7*t13*x4(t) - 92*x6(t)*t8 - x6(t)*t7*x4(t) + 276*x10(t)*x1(t)*t8*t13^2 + 276*x10(t)*t8*t13 + 276*x1(t)*t8*t13 + 276*t8)//(x10(t)*x1(t)*t13^2 + x10(t)*t13 + x1(t)*t13 + 1),
	x7'(t) = -92*x7(t)*t10 + x7(t)*x6(t)*t9 - 3*x7(t)*t9 + 15180*t10,
	x8'(t) = -x7(t)*t11 + 165*t11,
	x9'(t) = -2*x9(t)*u(t)*t12,
	x10'(t) = (-x8(t)*t16*x10(t) + x8(t)*t14 - t16*x10(t)*t15)//(x8(t) + t15),
	y1(t) = x3(t) + x1(t) + x4(t),
	y2(t) = -x9(t)*t18 + x5(t)*t18 + t18*x3(t) + t18*x4(t) + 1//3*t18,
	y3(t) = t19*x5(t) + t19*x4(t),
	y4(t) = -t20*x6(t) + 3*t20,
	y5(t) = x8(t)*t21,
	y6(t) = (x8(t)*t22*t17)//t11,
	y7(t) = x10(t),
	y8(t) = -x7(t) + 165
)
