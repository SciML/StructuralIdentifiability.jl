# TumorHu2019
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	z'(t) = (r3*k3*w(t) + r3*k3*m3 + r3*w(t)*v(t) + r3*m3*v(t) - k3*w(t)*mu3*z(t) - k3*mu3*m3*z(t) - w(t)*mu3*v(t)*z(t) - mu3*m3*v(t)*z(t) + v(t)*z(t)*b3)//(k3*w(t) + k3*m3 + w(t)*v(t) + m3*v(t)),
	v'(t) = (-k5*mu5*v(t) + b5*x(t)*z(t) - mu5*x(t)*v(t))//(k5 + x(t)),
	x'(t) = (-m1*r1*b1*x(t)^2 + m1*r1*x(t) - m1*b1^2*x(t)^2*y(t) + m1*b1*x(t)*y(t) - d1*x(t)*z(t) - w(t)*r1*b1*x(t)^2 + w(t)*r1*x(t) - w(t)*b1^2*x(t)^2*y(t) + w(t)*b1*x(t)*y(t))//(m1 + w(t)),
	w'(t) = (-m4*w(t)*mu4*k4 - m4*w(t)*mu4*x(t) + m4*x(t)*z(t)*b4 - w(t)*mu4*k4*v(t) - w(t)*mu4*x(t)*v(t) + k4*x(t)*r4*y(t) + x(t)^2*r4*y(t) + x(t)*v(t)*z(t)*b4)//(m4*k4 + m4*x(t) + k4*v(t) + x(t)*v(t)),
	y'(t) = (-mu2*w(t)*y(t) - mu2*k2*y(t) - b2^2*w(t)*y(t)^2 - b2*r2*w(t)*y(t)^2 - b2*r2*k2*y(t)^2 + b2*w(t)*y(t) + r2*w(t)*y(t) + r2*k2*y(t))//(w(t) + k2),
	y1(t) = z(t)
)
