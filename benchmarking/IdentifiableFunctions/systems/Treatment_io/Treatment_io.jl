# Treatment_io
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	S'(t) = (-b*In(t)*S(t) - b*S(t)*d*Tr(t))//N(t),
	In'(t) = (b*In(t)*S(t) + b*S(t)*d*Tr(t) - In(t)*N(t)*g - In(t)*N(t)*a)//N(t),
	Tr'(t) = In(t)*g - nu*Tr(t),
	N'(t) = 0,
	y1(t) = Tr(t),
	y2(t) = N(t)
)
