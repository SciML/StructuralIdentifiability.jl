z(t)diff(z(t), t) = c*w(t)*q*y(t) - h*z(t),
diff(v(t), t) = k*y(t) - v(t)*u,
diff(w(t), t) = -b*w(t) + c*w(t)*x(t)*y(t) - c*w(t)*q*y(t),
diff(y(t), t) = x(t)*v(t)*beta - a*y(t),
diff(x(t), t) = lm - x(t)*d - x(t)*v(t)*beta,
y1(t) = w(t),
y2(t) = z(t)
];
IdentifiabilityODE(sigma, GetParameters(sigma));