function model = SMTH()
syms v x z w y
syms b c h lm d k u q a beta
syms v0 x0 z0 w0 y0
model.sym.p = [b; c; h; lm; d; k; u; q; a; beta; v0; x0; z0; w0; y0];
model.sym.x = [v; x; z; w; y];
model.sym.g = [];
model.sym.x0 = [v0; x0; z0; w0; y0];
model.sym.xdot = [k*y - v*u
lm - x*d - x*v*beta
c*w*q*y - h*z
-b*w + c*w*x*y - c*w*q*y
x*v*beta - a*y];
model.sym.y = [z
w];
end