function model = SMTH()
syms S I W R
syms bi gam mu bw k xi a
syms S0 I0 W0 R0
model.sym.p = [bi; gam; mu; bw; k; xi; a; S0; I0; W0; R0];
model.sym.x = [S; I; W; R];
model.sym.g = [];
model.sym.x0 = [S0; I0; W0; R0];
model.sym.xdot = [-bi*S*I - S*mu - S*bw*W + R*a + mu
bi*S*I - gam*I + S*bw*W - mu*I
xi*I - xi*W
gam*I - R*mu - R*a];
model.sym.y = [k*I
S + R + I];
end