read "../IdentifiabilityODE.mpl";

sys := [
diff(v(t), t) = k*y(t) - v(t)*u,
diff(x(t), t) = lm - x(t)*d - x(t)*v(t)*beta,
diff(z(t), t) = c*w(t)*q*y(t) - h*z(t),
diff(w(t), t) = -b*w(t) + c*w(t)*x(t)*y(t) - c*w(t)*q*y(t),
diff(y(t), t) = x(t)*v(t)*beta - a*y(t),
y2(t) = z(t),
y1(t) = w(t)
];
CodeTools[CPUInfo](IdentifiabilityODE(sys, GetParameters(sys)));