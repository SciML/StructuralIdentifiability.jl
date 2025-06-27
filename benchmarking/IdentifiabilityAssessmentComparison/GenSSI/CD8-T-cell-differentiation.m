function model = SMTH()
syms N S M P E
syms rho_P mu_EE delta_NE mu_LE delta_EL mu_N mu_M delta_LM mu_PE mu_PL mu_LL mu_P rho_E
syms N0 S0 M0 P0 E0
model.sym.p = [rho_P; mu_EE; delta_NE; mu_LE; delta_EL; mu_N; mu_M; delta_LM; mu_PE; mu_PL; mu_LL; mu_P; rho_E; N0; S0; M0; P0; E0];
model.sym.x = [N; S; M; P; E];
model.sym.g = [];
model.sym.x0 = [N0; S0; M0; P0; E0];
model.sym.xdot = [-delta_NE*N*P - mu_N*N
-mu_LE*S*E + delta_EL*S - S^2*mu_LL - S*delta_LM
S*delta_LM - M*mu_M
rho_P*P^2 - S*P*mu_PL - E*mu_PE*P - P*mu_P
-mu_EE*E^2 + delta_NE*N*P - delta_EL*E + E*P*rho_E];
model.sym.y = [N
M
S + E];
end