r(t)diff(s(t), t) = -b1*b0*s(t)*x1(t)*i(t) - b0*s(t)*i(t) - s(t)*mu + mu + g*r(t),
diff(x2(t), t) = M*x1(t),
diff(i(t), t) = -nu*i(t) + b1*b0*s(t)*x1(t)*i(t) + b0*s(t)*i(t) - mu*i(t),
diff(r(t), t) = nu*i(t) - mu*r(t) - g*r(t),
diff(x1(t), t) = -M*x2(t),
y1(t) = i(t),
y2(t) = r(t)
];
IdentifiabilityODE(sigma, GetParameters(sigma));