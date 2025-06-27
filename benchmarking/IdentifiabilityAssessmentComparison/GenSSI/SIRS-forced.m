function model = SMTH()
syms s x2 i r x1
syms nu b1 b0 M mu g
syms s0 x20 i0 r0 x10
model.sym.p = [nu; b1; b0; M; mu; g; s0; x20; i0; r0; x10];
model.sym.x = [s; x2; i; r; x1];
model.sym.g = [];
model.sym.x0 = [s0; x20; i0; r0; x10];
model.sym.xdot = [-b1*b0*s*x1*i - b0*s*i - s*mu + mu + g*r
M*x1
-nu*i + b1*b0*s*x1*i + b0*s*i - mu*i
nu*i - mu*r - g*r
-M*x2];
model.sym.y = [i
r];
end