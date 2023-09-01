# cLV1 (1o)
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	pi2'(t) = u1(t)*B21*pi2(t)^2 + u1(t)*B21*pi2(t) - u1(t)*pi2(t)*pi1(t)*B11 + A23*pi2(t)^2*pi3(t) + A23*pi2(t)*pi3(t) + pi2(t)^3*A22 + pi2(t)^2*A21*pi1(t) + pi2(t)^2*g2 - pi2(t)^2*pi1(t)*A12 + pi2(t)^2*A22 + pi2(t)*A21*pi1(t) - pi2(t)*A11*pi1(t)^2 + pi2(t)*g2 - pi2(t)*pi1(t)*A13*pi3(t) - pi2(t)*pi1(t)*g1,
	pi3'(t) = u1(t)*B21*pi2(t)*pi3(t) + u1(t)*B31*pi3(t) - u1(t)*pi1(t)*pi3(t)*B11 + A23*pi2(t)*pi3(t)^2 + A33*pi3(t)^2 + pi2(t)^2*pi3(t)*A22 + pi2(t)*A21*pi1(t)*pi3(t) + pi2(t)*g2*pi3(t) - pi2(t)*pi1(t)*pi3(t)*A12 + pi2(t)*pi3(t)*A32 - A11*pi1(t)^2*pi3(t) + g3*pi3(t) + pi1(t)*A31*pi3(t) - pi1(t)*A13*pi3(t)^2 - pi1(t)*g1*pi3(t),
	pi1'(t) = u1(t)*B21*pi2(t)*pi1(t) - u1(t)*pi1(t)^2*B11 + u1(t)*pi1(t)*B11 + A23*pi2(t)*pi1(t)*pi3(t) + pi2(t)^2*pi1(t)*A22 + pi2(t)*A21*pi1(t)^2 + pi2(t)*g2*pi1(t) - pi2(t)*pi1(t)^2*A12 + pi2(t)*pi1(t)*A12 - A11*pi1(t)^3 + A11*pi1(t)^2 - pi1(t)^2*A13*pi3(t) - pi1(t)^2*g1 + pi1(t)*A13*pi3(t) + pi1(t)*g1,
	y1(t) = pi1(t)
)
