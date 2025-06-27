function model = SMTH()
syms S A E C Ninv J I
syms b alpha g2 k g1 q r
syms S0 A0 E0 C0 Ninv0 J0 I0
model.sym.p = [b; alpha; g2; k; g1; q; r; S0; A0; E0; C0; Ninv0; J0; I0];
model.sym.x = [S; A; E; C; Ninv; J; I];
model.sym.g = [];
model.sym.x0 = [S0; A0; E0; C0; Ninv0; J0; I0];
model.sym.xdot = [-b*S*Ninv*A*q - b*S*Ninv*I - b*S*Ninv*J
-E*k*r + E*k - A*g1
b*S*Ninv*A*q + b*S*Ninv*I + b*S*Ninv*J - E*k
alpha*I
0
alpha*I - g2*J
-alpha*I + E*k*r - g1*I];
model.sym.y = [Ninv
C];
end