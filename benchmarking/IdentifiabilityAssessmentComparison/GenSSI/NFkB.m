function model = SMTH()
syms x15 x13 x7 x9 x2 x8 x14 x5 x4 x12 x6 x1 x10 x11 x3
syms e2a i1 i1a c5 k_deg k2 t2 c3a k3 k1 c4a k_prod t1
syms x150 x130 x70 x90 x20 x80 x140 x50 x40 x120 x60 x10 x100 x110 x30
syms u
model.sym.p = [e2a; i1; i1a; c5; k_deg; k2; t2; c3a; k3; k1; c4a; k_prod; t1; x150; x130; x70; x90; x20; x80; x140; x50; x40; x120; x60; x10; x100; x110; x30];
model.sym.x = [x15; x13; x7; x9; x2; x8; x14; x5; x4; x12; x6; x1; x10; x11; x3];
model.sym.g = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; -k1 * x1; 0; 0; 0];
model.sym.x0 = [x150; x130; x70; x90; x20; x80; x140; x50; x40; x120; x60; x10; x100; x110; x30];
model.sym.xdot = [-1/2500*x15 + 1/2000000*x7
e2a*x14 + 1/2*x10*x6 - x2*x13 - 1/50000*x13
5*i1*x6 - 1/2*x11*x7
-1/2500*x9 - 1/2000000*x7
-1/5*x10*x2 - x2*x13 - x2*k_deg - x2*k2*u*x8 - x2*k3 + x1*u*k1 + x4*t1 + x5*t2
1/2*x9 - c5*x8
-5*e2a*x14 + 1/2*x11*x7
x2*x13 - x5*t2
1/5*x10*x2 - x4*t1
1/2000000*x7 - c3a*x12
-1/2*x10*x6 + 1/50000*x13 - i1*x6 - x5*t2
-x1*k_deg - x1*u*k1 + k_prod
-1/5*x10*x2 - x10*i1a - 1/2*x10*x6 - 1/10000*x10 + 1/2000*x11 + x12*c4a
5*x10*i1a - 1/2*x11*x7 - 1/400*x11
x2*k2*u*x8 + x2*k3 - x3*k_deg];
model.sym.y = [x2 + x3 + x1
x12
x10 + x13
x7
x9
x2];
end
