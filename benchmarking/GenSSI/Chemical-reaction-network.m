function model = SMTH()
syms x5 x6 x4 x2 x1 x3
syms k5 k3 k4 k2 k6 k1
syms x50 x60 x40 x20 x10 x30
model.sym.p = [k5; k3; k4; k2; k6; k1; x50; x60; x40; x20; x10; x30];
model.sym.x = [x5; x6; x4; x2; x1; x3];
model.sym.g = [];
model.sym.x0 = [x50; x60; x40; x20; x10; x30];
model.sym.xdot = [k5*x6 + k4*x6 - k6*x5*x3
-k5*x6 - k4*x6 + k6*x5*x3
-k3*x4 - k2*x4 + k1*x1*x2
k3*x4 + k2*x4 + k1*x1*x2
k4*x6 + k2*x4 - k1*x1*x2
k5*x6 + k3*x4 - k6*x5*x3];
model.sym.y = [x3
x2];
end