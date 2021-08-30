function model = SMTH()
syms x2 x4 x1 x3
syms b alpha c gama delta sigma beta
syms x20 x40 x10 x30
model.sym.p = [b; alpha; c; gama; delta; sigma; beta; x20; x40; x10; x30];
model.sym.x = [x2; x4; x1; x3];
model.sym.g = [];
model.sym.x0 = [x20; x40; x10; x30];
model.sym.xdot = [alpha*x1 - beta*x2
(gama*sigma*x2*x4 - delta*sigma*x3*x4) / (x3)
(-b*c*x1 - b*x1*x4 + 1) / (c + x4)
gama*x2 - delta*x3];
model.sym.y = [x1];
end