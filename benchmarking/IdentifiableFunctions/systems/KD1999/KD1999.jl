# KD1999
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	Cb'(t) = (-u1(t)*Cb(t) + V*Ca(t)*k0*Arr(t))//V,
	Tj'(t) = (Th*cph*u2(t)*roh - cph*u2(t)*roh*Tj(t) - UA*Tj(t) + UA*T(t))//(cph*Vh*roh),
	T'(t) = (-u1(t)*ro*T(t)*cp + u1(t)*ro*Ta*cp - V*Ca(t)*k0*Arr(t)*DH - UA*Tj(t) + UA*T(t))//(V*ro*cp),
	Ca'(t) = (-u1(t)*Ca(t) + u1(t)*Ca0 - V*Ca(t)*k0*Arr(t))//V,
	Arr'(t) = (-u1(t)*ro*E*Arr(t)*T(t)*cp + u1(t)*ro*E*Arr(t)*Ta*cp - V*Ca(t)*k0*E*Arr(t)^2*DH - UA*Tj(t)*E*Arr(t) + UA*E*Arr(t)*T(t))//(V*ro*R*T(t)^2*cp),
	y2(t) = T(t),
	y1(t) = Cb(t)
)
