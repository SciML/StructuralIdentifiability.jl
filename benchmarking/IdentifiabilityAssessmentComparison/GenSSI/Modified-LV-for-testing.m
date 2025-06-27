function model = SMTH()
syms x1 x2
syms a b d c
syms x10 x20
model.sym.p = [a; b; d; c; x10; x20];
model.sym.x = [x1; x2];
model.sym.g = [];
model.sym.x0 = [x10; x20];
model.sym.xdot = [a*x1 + b*x1 - x2*c*x1
-a*b*x2 + d*x2*x1];
model.sym.y = [x1];
end