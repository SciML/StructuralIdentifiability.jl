## Benchmark results

2026-07-28T01:16:05.923

- Benchmarked function: `find_identifiable_functions`
- Workers: 30
- Timeout: 7200 s
- Memory limit: 28192 MiB

|Model|funcs (cmp_default)|funcs (cmp_lie)|
|---|---|---|
|Akt pathway|[reaction_8_k1, reaction_7_k1, reaction_6_k1, reaction_5_k2, reaction_4_k1, reaction_3_k1, reaction_2_k2, reaction_1_k1 - reaction_1_k2 - reaction_9_k1, reaction_2_k1//reaction_5_k1, a3//reaction_5_k1, a2//reaction_5_k1, a1//reaction_5_k1, pS6(t)*reaction_5_k1, pAkt_S6(t)*reaction_5_k1, S6(t)*reaction_5_k1, pAkt(t)*reaction_5_k1, Akt(t)*reaction_5_k1, pEGFR_Akt(t)*reaction_5_k1, pEGFR(t)*reaction_5_k1, EGF_EGFR(t)*reaction_5_k1*reaction_9_k1]|[reaction_8_k1, reaction_7_k1, reaction_6_k1, reaction_5_k2, reaction_4_k1, reaction_3_k1, reaction_2_k2, reaction_1_k1 - reaction_1_k2 - reaction_9_k1, reaction_2_k1//reaction_5_k1, a3//reaction_5_k1, a2//reaction_5_k1, a1//reaction_5_k1, pS6(t)*reaction_5_k1, pAkt(t)*a2 + pAkt_S6(t)*a2, pEGFR(t)*a1 + pEGFR_Akt(t)*a1, EGF_EGFR(t)*reaction_5_k1*reaction_9_k1, S6(t)*reaction_5_k1, pAkt_S6(t)*reaction_5_k1, Akt(t)*reaction_5_k1, pEGFR_Akt(t)*reaction_5_k1]|
|Bilirubin2_io|[k01, k21 + k31 + k41, k12 + k13 + k14, k21*k31 + k21*k41 + k31*k41, k12*k13 + k12*k14 + k13*k14, k12*k31 + k12*k41 + k13*k21 + k13*k41 + k14*k21 + k14*k31, k21*k31*k41, k12*k13*k14, x1(t), x2(t) + x3(t) + x4(t), x2(t)*x3(t) + x2(t)*x4(t) + x3(t)*x4(t), x2(t)*k31 + x2(t)*k41 + x3(t)*k21 + x3(t)*k41 + x4(t)*k21 + x4(t)*k31, x2(t)*k13 + x2(t)*k14 + x3(t)*k12 + x3(t)*k14 + x4(t)*k12 + x4(t)*k13]|[k01, k21 + k31 + k41, k12 + k13 + k14, k21*k31 + k21*k41 + k31*k41, k12*k13 + k12*k14 + k13*k14, k12*k31 + k12*k41 + k13*k21 + k13*k41 + k14*k21 + k14*k31, k21*k31*k41, k12*k13*k14, x2(t) + x3(t) + x4(t), x1(t), x2(t)*k31 + x2(t)*k41 + x3(t)*k21 + x3(t)*k41 + x4(t)*k21 + x4(t)*k31, x2(t)*k13 + x2(t)*k14 + x3(t)*k12 + x3(t)*k14 + x4(t)*k12 + x4(t)*k13]|
|Biohydrogenation_io|[k7, k6, k5, k10 + 2*k8, k9^2, k10*k9, x5(t), x4(t), x6(t) + k8]|[k7, k6, k5, k10 + 2*k8, k9^2, k10*k9, x4(t), x5(t), x6(t) - 1//2*k10]|
|Bruno2016|[kbeta10, kbeta, kcryOH + kcrybeta, beta10(t), beta(t), cry(t)*kcryOH]|[kbeta10, kbeta, kcryOH + kcrybeta, beta(t), beta10(t), cry(t)*kcryOH]|
|CD8 T cell differentiation|[mu_PL, mu_PE, mu_P, mu_N, mu_M, mu_LL, mu_LE, mu_EE, delta_LM, delta_EL, rho_E//rho_P, delta_NE//rho_P, M(t), S(t), E(t), N(t), P(t)*rho_P]|[mu_PL, mu_PE, mu_P, mu_N, mu_M, mu_LL, mu_LE, mu_EE, delta_LM, delta_EL, rho_E//rho_P, delta_NE//rho_P, M(t), N(t), S(t), E(t), P(t)*rho_P]|
|CGV1990|[k7, k6, k4, k3, R*V36 + 1//25*S*V3, R*V3 + S*V36, R*V36*k5 + 1//5*S*V36*k5, (V3*k5 + 5*V36*k5)//V3, q7(t), q3(t), q1(t), q35(t) + q36(t), q35(t)*q36(t), q35(t)*S*V36 + q36(t)*R*V3]|[k7, k6, k4, k3, R*V36 + 1//25*S*V3, R*V3 + S*V36, R*V36*k5 + 1//5*S*V36*k5, (V3*k5 + 5*V36*k5)//V3, q7(t), q1(t), q35(t)*S*V36 + q36(t)*R*V3, q35(t)*R*V36 + 1//25*q36(t)*S*V3, q3(t)]|
|Chemical reaction network|[k6, k5, k4, k3, k2, k1, x6(t), x5(t), x4(t), x3(t), x2(t), x1(t)]|[k6, k5, k4, k3, k2, k1, x6(t), x5(t), x3(t), x1(t), x2(t), x4(t)]|
|Crauste_SI|[mu_PL, mu_PE, mu_P, mu_N, mu_M, mu_LL, mu_LE, mu_EE, delta_LM, delta_EL, rho_E//rho_P, delta_NE//rho_P, M(t), S(t), E(t), N(t), P(t)*rho_P]|[mu_PL, mu_PE, mu_P, mu_N, mu_M, mu_LL, mu_LE, mu_EE, delta_LM, delta_EL, rho_E//rho_P, delta_NE//rho_P, M(t), N(t), S(t), E(t), P(t)*rho_P]|
|Fujita|[reaction_8_k1, reaction_7_k1, reaction_6_k1, reaction_5_k2, reaction_4_k1, reaction_3_k1, reaction_2_k2, reaction_1_k1 - reaction_1_k2 - reaction_9_k1, reaction_2_k1//reaction_5_k1, a3//reaction_5_k1, a2//reaction_5_k1, a1//reaction_5_k1, pS6(t)*reaction_5_k1, pAkt_S6(t)*reaction_5_k1, S6(t)*reaction_5_k1, pAkt(t)*reaction_5_k1, Akt(t)*reaction_5_k1, pEGFR_Akt(t)*reaction_5_k1, pEGFR(t)*reaction_5_k1, EGF_EGFR(t)*reaction_5_k1*reaction_9_k1]|[reaction_8_k1, reaction_7_k1, reaction_6_k1, reaction_5_k2, reaction_4_k1, reaction_3_k1, reaction_2_k2, reaction_1_k1 - reaction_1_k2 - reaction_9_k1, reaction_2_k1//reaction_5_k1, a3//reaction_5_k1, a2//reaction_5_k1, a1//reaction_5_k1, pS6(t)*reaction_5_k1, pAkt(t)*a2 + pAkt_S6(t)*a2, pEGFR(t)*a1 + pEGFR_Akt(t)*a1, EGF_EGFR(t)*reaction_5_k1*reaction_9_k1, S6(t)*reaction_5_k1, pAkt_S6(t)*reaction_5_k1, Akt(t)*reaction_5_k1, pEGFR_Akt(t)*reaction_5_k1]|
|Goodwin oscillator|[sigma, c, b, beta + delta, beta*delta, x4(t), x1(t), (alpha*gama)//x3(t), (x2(t)*gama - x3(t)*delta)//x3(t)]|[sigma, c, b, beta + delta, beta*delta, x1(t), x4(t), (x2(t)*gama + 1//2*x3(t)*beta - 1//2*x3(t)*delta)//(alpha*gama), (alpha*gama)//x3(t)]|
|HIV|[u, h, d, b, a, lm//q, c*q^2, beta*k*q, z(t), w(t), v(t)*beta, y(t)//q, x(t)//q]|[u, h, d, b, a, lm//q, c*q^2, beta*k*q, v(t)*beta, z(t), y(t)//q, w(t), x(t)//q]|
|HIV2_io|[s, d, b, c + k1 + w1 + w2, k2*q2, c*k1 + c*w1 + c*w2 + k1*w2 + w1*w2, c*k1*w2 + c*w1*w2, k1*k2*q1 + k1*k2*q2 + k2*q2*w1, x4(t), x1(t), x3(t)*k2 - x4(t)*c, x2(t)*k1*k2 + x3(t)*k1*k2 + x3(t)*k2*w1 + x4(t)*k1*w2 + x4(t)*w1*w2]|[s, d, b, c + k1 + w1 + w2, k2*q2, c*k1 + c*w1 + c*w2 + k1*w2 + w1*w2, c*k1*w2 + c*w1*w2, k1*k2*q1 + k1*k2*q2 + k2*q2*w1, x4(t), x1(t), x3(t)*k2 - x4(t)*c, x2(t)*k1*k2 + x3(t)*k1*k2 + x3(t)*k2*w1 + x4(t)*k1*w2 + x4(t)*w1*w2]|
|HighDimNonLin|[vm, p9, p8, p7, p6, p5, p4, p3, p20, p2, p19, p18, p17, p16, p15, p14, p13, p12, p11, p10, p1, km, x20(t), x19(t), x18(t), x17(t), x16(t), x15(t), x14(t), x13(t), x12(t), x11(t), x10(t), x9(t), x8(t), x7(t), x6(t), x5(t), x4(t), x3(t), x2(t), x1(t)]|[vm, p9, p8, p7, p6, p5, p4, p3, p20, p2, p19, p18, p17, p16, p15, p14, p13, p12, p11, p10, p1, km, x20(t), x19(t), x18(t), x17(t), x16(t), x15(t), x14(t), x13(t), x12(t), x11(t), x10(t), x9(t), x8(t), x7(t), x6(t), x5(t), x4(t), x3(t), x2(t), x1(t)]|
|JAK-STAT 1|[t9, t8, t7, t6, t5, t4, t3, t20, t2, t19, t18, t16, t14, t13, t12, t10, t1, t17*t22, t15*t21, t11*t21, x10(t), x9(t), x7(t), x6(t), x5(t), x4(t), x3(t), x2(t), x1(t), x8(t)*t21]|[t9, t8, t7, t6, t5, t4, t3, t20, t2, t19, t18, t16, t14, t13, t12, t10, t1, t17*t22, t15*t21, t11*t21, x5(t), x2(t), x1(t) + x3(t) + x4(t), x9(t), (x8(t)*t17*t22)//t11, x4(t), x3(t), x7(t), x10(t), x6(t)]|
|KD1999|[Vh, V, Th, Ta, Ca0, DH//UA, E//R, (cph*roh)//UA, (cp*ro)//UA, Tj(t), T(t), Cb(t), Ca(t), Arr(t)*k0]|[Vh, V, Th, Ta, Ca0, DH//UA, E//R, (cph*roh)//UA, (cp*ro)//UA, Cb(t), Ca(t), Tj(t), T(t), Arr(t)*k0]|
|LLW1987_io|[p1 + p3, p2*p4, p1*p3, x3(t), x1(t)*x2(t), x1(t)*p4 + x2(t)*p2, x1(t)*p1*p4 + x2(t)*p2*p3]|[p1 + p3, p2*p4, p1*p3, x1(t)*p4 + x2(t)*p2, x3(t), x1(t)*x2(t), (x1(t)*p4 - x2(t)*p2)//(p1 - p3)]|
|LeukaemiaLeon2021|-|-|
|MAPK model (5 outputs bis)|-|-|
|MAPK model (5 outputs)|[gamma1110, gamma1101, gamma1100, gamma1000, gamma0100, c1011, c0111, c0011, c0010, c0001, beta11, beta10, beta01, b10, b01, b00, alpha11, alpha10, alpha01, a10, a01, a00, S11(t), S10(t), S01(t), S00(t), F(t), K(t), FS11(t), FS10(t), FS01(t), KS10(t), KS01(t), KS00(t)]|[gamma1110, gamma1101, gamma1100, gamma1000, gamma0100, c1011, c0111, c0011, c0010, c0001, beta11, beta10, beta01, b10, b01, b00, alpha11, alpha10, alpha01, a10, a01, a00, S00(t), F(t), K(t), FS11(t), S01(t), FS01(t), KS00(t), FS10(t), KS10(t), KS01(t), S11(t), S10(t)]|
|MAPK model (6 outputs)|[gamma1110, gamma1101, gamma1100, gamma1000, gamma0100, c1011, c0111, c0011, c0010, c0001, beta11, beta10, beta01, b10, b01, b00, alpha11, alpha10, alpha01, a10, a01, a00, S11(t), S10(t), S01(t), S00(t), F(t), K(t), FS11(t), FS10(t), FS01(t), KS10(t), KS01(t), KS00(t)]|[gamma1110, gamma1101, gamma1100, gamma1000, gamma0100, c1011, c0111, c0011, c0010, c0001, beta11, beta10, beta01, b10, b01, b00, alpha11, alpha10, alpha01, a10, a01, a00, S00(t), F(t), K(t), FS11(t), S01(t), FS01(t), KS00(t), FS10(t), KS10(t), KS01(t), S11(t), S10(t)]|
|Modified LV for testing|[d, a + b, a*b, x1(t), x2(t)*c]|[d, a + b, a*b, x1(t), x2(t)*c]|
|NFkB|-|-|
|PK1|[k6, k5, k4, k3 + k7, k1 + k2, k2*s3, k1*s2, k1*k3*s3, x1(t), x4(t)*s2, x3(t)*s3, x2(t)*s2]|[k6, k5, k4, k3 + k7, k1 + k2, k2*s3, k1*s2, k1*k3*s3, x1(t), x4(t)*s2, x3(t)*s3, x2(t)*s2]|
|PK2|[n, kc, ka, b2, b1, a2, a1, x3(t), x2(t), x1(t), x0(t)]|[n, kc, ka, b2, b1, a2, a1, x3(t), x1(t), x2(t), x0(t)]|
|Pharm|[n, kc, ka, b2, b1, a2, a1, x3(t), x2(t), x1(t), x0(t)]|[n, kc, ka, b2, b1, a2, a1, x3(t), x1(t), x2(t), x0(t)]|
|Phosphorylation|[k6, k5, k4, k3, k2, k1, x3(t), x1(t), x2(t), x4(t), x6(t), x5(t)]|[k6, k5, k4, k3, k2, k1, x4(t), x2(t), x1(t), x3(t), x6(t), x5(t)]|
|Pivastatin|[r3, r1, k4, k3, k2, T0*k1, T0*k, x3(t)*k1, x2(t)*k1, x1(t)*k1]|[r3, r1, k4, k3, k2, T0*k1, T0*k, x3(t)*k1, x2(t)*k1, x1(t)*k1]|
|QWWC|-|-|
|QY|[alpa, Mar, Ks, M + siga1 + siga2, phi//(Mar - siga1), M*siga1 + M*siga2 + siga1*siga2, M*siga1*siga2, (M*phi*siga2 - M*siga2)//beta_SA, (Ks*beta_SA + M*beta_SA - 1//3*Mar*beta_SA*phi + 1//3*Mar*beta_SA + 1//3*beta_SA*phi*siga2 + 2//3*beta_SA*siga1 + 2//3*beta_SA*siga2 - 1//3*beta_SI*phi*siga2 + 1//3*beta_SI*siga2)//beta_SA, P4(t), P3(t), P2(t), P1(t), P0(t), P5(t)*Mar + beta]|[alpa, Mar, Ks, M + siga1 + siga2, phi//(Mar - siga1), M*siga1 + M*siga2 + siga1*siga2, M*siga1*siga2, (M*phi*siga2 - M*siga2)//beta_SA, (Ks*beta_SA + M*beta_SA - 1//3*Mar*beta_SA*phi + 1//3*Mar*beta_SA + 1//3*beta_SA*phi*siga2 + 2//3*beta_SA*siga1 + 2//3*beta_SA*siga2 - 1//3*beta_SI*phi*siga2 + 1//3*beta_SI*siga2)//beta_SA, P3(t), P2(t), P1(t), P0(t), (P5(t)*Mar + beta)//Mar, P4(t)]|
|Ruminal lipolysis|[k4, k3, k2, x5(t), x4(t), x3(t), x2(t), x1(t)]|[k4, k3, k2, x5(t), x3(t), x1(t), x2(t), x4(t)]|
|SEAIJRC Covid model|[k, g2, g1, b, alpha, (q*r - q)//r, Ninv(t), C(t), J(t), I(t), A(t)*q, E(t)*r, S(t)*r]|[k, g2, g1, b, alpha, Ninv(t), (q*r - q)//r, C(t), J(t), I(t), A(t)*q, S(t)*r, E(t)*r]|
|SEIR 34|[mu, epsilon + gamma, epsilon*gamma, K*epsilon, beta*epsilon*r, S(t), N(t), A(t), E(t) + I(t), I(t)*gamma]|[mu, N(t), A(t), epsilon + gamma, epsilon*gamma, K*epsilon, beta*epsilon*r, I(t)*gamma, S(t), E(t) + I(t)]|
|SEIR 36 ref|[s_d, s, phi_e, phi, mu_i, mu_d, mu_0, gamma_d, gamma, beta_d, beta, F(t), Di(t), De(t), I(t), E(t), S(t), q(t), nu(t), N(t)]|[s_d, s, phi_e, phi, mu_i, mu_d, mu_0, gamma_d, gamma, beta_d, beta, q(t), nu(t), N(t), F(t), De(t), Di(t), I(t), S(t), E(t)]|
|SEIR2T|[nu, b, a, Cu(t), N(t), In(t), E(t), S(t)]|[nu, b, a, N(t), Cu(t), In(t), S(t), E(t)]|
|SEIRT|[beta, alpha + lambda, alpha*lambda, N(t), I(t), S(t)*alpha, E(t)*alpha - I(t)*lambda]|[beta, N(t), alpha + lambda, alpha*lambda, I(t), S(t)*alpha, (E(t) + I(t))//lambda]|
|SEIR_1_io|[gamma, beta//psi, gamma*psi - psi*v, psi*v - psi - v, Q(t), I(t)*psi, S(t)*psi*v, E(t)*psi*v + I(t)*psi*v]|[gamma, beta//psi, gamma*psi - psi*v, psi*v - psi - v, Q(t), I(t)*psi, S(t)*psi*v, E(t)*psi*v + I(t)*psi*v]|
|SEUIR|[z, d, (N*w)//beta, I(t), E(t)*w, S(t)*w, U(t)*w + I(t)*w]|[z, d, (N*w)//beta, I(t), U(t)*w + I(t)*w, S(t)*w, E(t)*w]|
|SIR 19|[r, q, pp, mu, beta, D(t), C(t), I(t), S(t), N(t)]|[r, q, pp, mu, beta, N(t), D(t), C(t), I(t), S(t)]|
|SIR 21|[r, q, pp, mu, beta, D(t), C(t), I(t), S(t), N(t)]|[r, q, pp, mu, beta, N(t), D(t), C(t), I(t), S(t)]|
|SIR 24|[K, gamma + mu, c*phi, I(t), S(t), A(t) - mu]|[K, gamma + mu, A(t) - mu, c*phi, I(t), S(t)]|
|SIR 6|[gamma, K//beta, N(t), I(t)*beta, S(t)*beta]|[gamma, N(t), K//beta, S(t)*beta, I(t)*beta]|
|SIRS forced|[nu, mu, g, b0, M^2, r(t), i(t), s(t), x1(t)*b1, x2(t)*M*b1]|[nu, mu, g, b0, M^2, (x2(t)*b1)//M, r(t), x1(t)*b1, i(t), s(t)]|
|SIWR original|[xi, mu, k, gam, bw, bi, a, R(t), W(t), I(t), S(t)]|[xi, mu, k, gam, bw, bi, a, W(t), R(t), I(t), S(t)]|
|SIWR with extra output|[xi, mu, k, gam, bw, bi, a, R(t), W(t), I(t), S(t)]|[xi, mu, k, gam, bw, bi, a, W(t), R(t), S(t) + I(t) + R(t), I(t)]|
|SLIQR|[s, b, Ninv, a + g, a*e*g, a*g + e*g*s - g*s, In(t), S(t)*a, Q(t)*a - Q(t)*s, L(t)*a - In(t)*g + Q(t)*s]|[s, b, Ninv, a + g, a*e*g, a*g + e*g*s - g*s, In(t), Q(t)*a - Q(t)*s, S(t)*a, L(t)*a - In(t)*g + Q(t)*s]|
|St|[T, Dd, T*a + T*d + T*dr + e + g - r - rR, (d*rR - dr*r)//(d - dr), (a*r - a*rR + d*e + d*g - dr*e - dr*g)//(d - dr), (a^2 + 2*a*d + d^2 + dr^2)//(a*dr + d*dr), (a*dr*r - a*dr*rR + d^2*g + d*dr*e - d*dr*g - dr^2*e)//(a*d - a*dr + d^2 - dr^2), S(t) + R(t), W(t)*a + W(t)*d + W(t)*dr + e + g - r - rR, S(t)*T*d - S(t)*r + R(t)*T*dr - R(t)*rR]|[T, Dd, T*a + T*d + T*dr + e + g - r - rR, (d*rR - dr*r)//(d - dr), (a*r - a*rR + d*e + d*g - dr*e - dr*g)//(d - dr), (a^2 + 2*a*d + d^2 + dr^2)//(a*dr + d*dr), (a*dr*r - a*dr*rR + d^2*g + d*dr*e - d*dr*g - dr^2*e)//(a*d - a*dr + d^2 - dr^2), S(t) + R(t), W(t)*a + W(t)*d + W(t)*dr + e + g - r - rR, (S(t)*d + R(t)*dr)//(a + d + dr)]|
|Transfection_4State|[d3, d1, b, d2//kTL, GFP(t), mRNAenz(t)*kTL, enz(t)*kTL, mRNA(t)*kTL]|[d3, d1, b, d2//kTL, GFP(t), mRNAenz(t)*kTL, enz(t)*kTL, mRNA(t)*kTL]|
|Treatment_io|[a + g + nu, b//g, d*g + nu, a*nu + g*nu, N(t), Tr(t), S(t)*g, In(t)*g - Tr(t)*nu]|[N(t), a + g + nu, b//g, d*g + nu, a*nu + g*nu, Tr(t), S(t)*g, In(t)*g - Tr(t)*nu]|
|TumorHu2019|-|-|
|TumorPillis2007|-|-|
|cLV1 (1o)|-|-|
|cLV1 (2o)|[g3, g2, g1, B31, B21, B11, A32, A31, A22, A21, A12, A11, A23//A33, A13//A33, pi2(t), pi1(t), pi3(t)*A33]|[g3, g2, g1, B31, B21, B11, A32, A31, A22, A21, A12, A11, A23//A33, A13//A33, pi2(t), pi1(t), pi3(t)*A33]|
|cholera|[xi, mu, k, gam, bw, bi, a, R(t), W(t), I(t), S(t)]|[xi, mu, k, gam, bw, bi, a, W(t), R(t), S(t) + I(t) + R(t), I(t)]|
|generalizedLoktaVolterra (1o)|[r2, r1, beta21, beta11, beta12//beta22, x1(t), x2(t)*beta22]|[r2, r1, beta21, beta11, beta12//beta22, x1(t), x2(t)*beta22]|
|generalizedLoktaVolterra (2o)|[r2, r1, beta22, beta21, beta12, beta11, x2(t), x1(t)]|[r2, r1, beta22, beta21, beta12, beta11, x2(t), x1(t)]|
|p53|[p9, p8, p7, p6, p5, p4, p3, p25, p24, p23, p21, p20, p18, p17, p16, p15, p14, p13, p12, p11, p10, p1, p22^4, x4(t), x3(t), x2(t), x1(t)]|[p9, p8, p7, p6, p5, p4, p3, p25, p24, p23, p21, p20, p18, p17, p16, p15, p14, p13, p12, p11, p10, p1, p22^4, x3(t), x2(t), x1(t), x4(t)]|

*Benchmarking environment:*

* Total RAM (GiB): 251
* Processor: INTEL(R) XEON(R) GOLD 6538Y+
* Julia version: 1.12.6

Versions of the dependencies:

* AbstractAlgebra : 0.50.1
* Combinatorics : 1.1.0
* DataStructures : 0.19.6
* Dates : 1.11.0
* Groebner : 0.10.4
* IterTools : 1.10.0
* LinearAlgebra : 1.12.0
* Logging : 1.11.0
* MacroTools : 0.5.16
* Nemo : 0.56.1
* ParamPunPam : 0.5.8
* PrecompileTools : 1.3.4
* Primes : 0.5.7
* Random : 1.11.0
* RationalFunctionFields : 0.3.2
* TimerOutputs : 0.5.29
