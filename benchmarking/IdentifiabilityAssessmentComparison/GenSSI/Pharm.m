function model = SMTH()
syms x3 x1 x2 x0
syms a2 ka n b2 kc b1 a1
syms x30 x10 x20 x00
model.sym.p = [a2; ka; n; b2; kc; b1; a1; x30; x10; x20; x00];
model.sym.x = [x3; x1; x2; x0];
model.sym.g = [];
model.sym.x0 = [x30; x10; x20; x00];
model.sym.xdot = [-b2*x3 + b2*x2
a2*x0 - a2*x1
(ka*kc*b1*x3 - ka*kc*b1*x2 + ka*b1*x0*x3 - ka*b1*x0*x2 - n*kc*x2 + kc*b1*x3*x2 - kc*b1*x2^2) / (ka*kc + ka*x0 + kc*x2)
(-ka*n*x0 - ka*kc*a1*x0 + ka*kc*a1*x1 - ka*a1*x0^2 + ka*a1*x0*x1 - kc*a1*x0*x2 + kc*a1*x1*x2) / (ka*kc + ka*x0 + kc*x2)];
model.sym.y = [x0];
end