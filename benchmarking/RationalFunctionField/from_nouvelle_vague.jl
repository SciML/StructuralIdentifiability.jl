# Examples arising from the wave-approach for TumorPillis

R, (T, L, N, C, I, M, KC, KL, KN, KT, a, alpha1, alpha2, b, beta, c1, f, g, gI, gamma, gt, h, m, muI, p, pI, pt, q, r2, ucte, w) = polynomial_ring(QQ, ["T", "L", "N", "C", "I", "M", "KC", "KL", "KN", "KT", "a", "alpha1", "alpha2", "b", "beta", "c1", "f", "g", "gI", "gamma", "gt", "h", "m", "muI", "p", "pI", "pt", "q", "r2", "ucte", "w"])
 
push!(
    cases_simplification,
    Dict(
        :description => "New tide for TumorPillis, step 1",
        :gens => [M, N, L, -M*gamma, -gamma, M*gamma^2, -N*KN // one(R)],
    )
)

push!(
    cases_simplification,
    Dict(
        :description => "New tide for TumorPillis, step 2",
        :gens => [gamma, KN, M, N, L, T*L*q - T*C*r2, T*L*q - T*C*r2, (L*gI*pI)//(I^2 + 2*I*gI + gI^2), (L*gI*pI)//(I^2 + 2*I*gI + gI^2)],
    )
)

push!(
    cases_simplification,
    Dict(
        :description => "New tide for TumorPillis, step 3",
        :gens => [gamma, KN, M, N, L, T*L*q - T*C*r2, (gI*pI)//(I^2 + 2*I*gI + gI^2), T*q, T*q, T*q, -T*q, -T*q, (-1//2*I - 1//2*gI)//L, (-1//2*I - 1//2*gI)//L, L*w - muI, (T*L*gI*pI*q)//(I^2 + 2*I*gI + gI^2), (-T*L*gt*pt)//(T^2 + 2*T*gt + gt^2)],
    )
)

push!(
    cases_simplification,
    Dict(
        :description => "New tide for TumorPillis, step 4",
        :gens => [gamma, KN, M, N, L, gI*pI, T*q, I + gI, L*w - muI, (C*r2)//q, (T*gt*pt)//(T^2 + 2*T*gt + gt^2), w, w, -T*KT*q, (-C*KC*r2)//q, T*L*q*w - T*C*r2*w, (-C*M*KC*r2 - C*beta*r2 + alpha2*r2)//q, (L*gI*pI*w)//(I^2 + 2*I*gI + gI^2), -T^2*a*b*q - T*N*c1*q - T*M*KT*q + T*a*q, -T^2*a*b*q - T*N*c1*q - T*M*KT*q + T*a*q, -T^2*a*b*q - T*N*c1*q - T*M*KT*q + T*a*q, -T^2*a*b*q - T*N*c1*q - T*M*KT*q + T*a*q, T^2*a*b*q + T*N*c1*q + T*M*KT*q - T*a*q, T^2*a*b*q + T*N*c1*q + T*M*KT*q - T*a*q, T^2*a*b*q + T*N*c1*q + T*M*KT*q - T*a*q, 2*T^2*a*b*q + T*N*c1*q + T*M*KT*q - T*a*q, -T^2*a*b*q - T*N*c1*q - T*M*KT*q + T*a*q],
    )
)

push!(
    cases_simplification,
    Dict(
        :description => "New tide for TumorPillis, step 5",
        :gens => [w, muI, gamma, KT, KN, KC, M, N, L, gI*pI, T*q, T*C*r2, I + gI, N*c1 - a, T*C*beta*r2 - T*alpha2*r2, (a*b)//(N*q), (T*gt*pt)//(T^2 + 2*T*gt + gt^2), -N*KN*c1, -T*C*KC*beta*r2 - T*C*KT*beta*r2 + T*KT*alpha2*r2, (C*M^2*KC^2*r2 + 2*C*M*KC*beta*r2 + C*M*KC*gamma*r2 + C*beta^2*r2 - M*KC*alpha2*r2 - alpha2*beta*r2)//q, (T*I*w + T*pt + I*gt*w)//(T + gt), (T*I*w + T*pt + I*gt*w)//(T + gt), (C*M^2*KC^2*r2 + 2*C*M*KC*beta*r2 + C*M*KC*gamma*r2 + C*beta^2*r2 - M*KC*alpha2*r2 - alpha2*beta*r2)//q, (T*I*w + T*pt + I*gt*w)//(T + gt), (T^2*gt*pt - T*gt^2*pt)//(T^3 + 3*T^2*gt + 3*T*gt^2 + gt^3), (T^2*KT*gt*pt - T*KT*gt^2*pt)//(T^3 + 3*T^2*gt + 3*T*gt^2 + gt^3), (T^2*gt*pt - T*gt^2*pt)//(T^3 + 3*T^2*gt + 3*T*gt^2 + gt^3), (-T^2*L*gt*pt + T*L*gt^2*pt)//(T^3 + 3*T^2*gt + 3*T*gt^2 + gt^3), (T^2*KT*gt*pt - T*KT*gt^2*pt)//(T^3 + 3*T^2*gt + 3*T*gt^2 + gt^3), (T*L*I*w + T*L*pt - T*I*muI + L*I*gt*w - I*gt*muI)//(T + gt), (T*L*I*w + T*L*pt - T*I*muI + L*I*gt*w - I*gt*muI)//(T + gt), (T*L*I*w + T*L*pt - T*I*muI + L*I*gt*w - I*gt*muI)//(T + gt), -T^2*C*a*b*beta*r2 + T^2*a*alpha2*b*r2 - T*N*C*beta*c1*r2 + T*N*alpha2*c1*r2 - T*C*M*KC*beta*r2 - T*C*M*KT*beta*r2 + T*C*a*beta*r2 - T*C*beta^2*r2 + T*M*KT*alpha2*r2 - T*a*alpha2*r2 + T*alpha2*beta*r2, T^2*C*a*b*beta*r2 - T^2*a*alpha2*b*r2 + T*N*C*beta*c1*r2 - T*N*alpha2*c1*r2 + T*C*M*KC*beta*r2 + T*C*M*KT*beta*r2 - T*C*a*beta*r2 + T*C*beta^2*r2 - T*M*KT*alpha2*r2 + T*a*alpha2*r2 - T*alpha2*beta*r2, 2*T^2*C*a*b*beta*r2 - 2*T^2*a*alpha2*b*r2 + T*N*C*beta*c1*r2 - T*N*alpha2*c1*r2 + T*C*M*KC*beta*r2 + T*C*M*KT*beta*r2 - T*C*a*beta*r2 + T*C*beta^2*r2 - T*M*KT*alpha2*r2 + T*a*alpha2*r2 - T*alpha2*beta*r2, (T^3*N*p + 2*T^2*N*h*p - T*N*g*h + T*N*h^2*p)//(T^2 + 2*T*h + h^2), (T^3*N*c1*p + 2*T^2*N*c1*h*p - T*N*c1*g*h + T*N*c1*h^2*p)//(T^2 + 2*T*h + h^2), (-T^3*a*b*p - 2*T^2*a*b*h*p + T*a*b*g*h - T*a*b*h^2*p)//(T^2*N*q + 2*T*N*h*q + N*h^2*q), (T^3*N*p + 2*T^2*N*h*p - T*N*g*h + T*N*h^2*p)//(T^2 + 2*T*h + h^2), (T^3*N*p + 2*T^2*N*h*p - T*N*g*h + T*N*h^2*p)//(T^2 + 2*T*h + h^2), (T^3*N*p + 2*T^2*N*h*p - T*N*g*h + T*N*h^2*p)//(T^2 + 2*T*h + h^2), (T^3*N*p + 2*T^2*N*h*p - T*N*g*h + T*N*h^2*p)//(T^2 + 2*T*h + h^2)],
    )
)

push!(
    cases_simplification,
    Dict(
        :description => "New tide for TumorPillis, step 6",
        :gens => [w, pt, pI, muI, gamma, gI, c1, beta, a, KT, KN, KC, M, I, N, L, gt*q, b*gt, T*q, C//alpha2, (alpha2*r2)//q, (T^2*p + 2*T*h*p - g*h + h^2*p)//(T^2*q + 2*T*h*q + h^2*q), (T^3*q + 3//2*T^2*h*q - 1//2*h^3*q)//(T*g*h), (-1//2*T^3*q - 3//2*T^2*h*q - 3//2*T*h^2*q - 1//2*h^3*q)//(T*g*h), (T^4*a*b*q + T^3*N*c1*q + T^3*M*KT*q + 3//2*T^3*a*b*h*q - T^3*a*q + 3//2*T^2*N*c1*h*q + 3//2*T^2*M*KT*h*q - 3//2*T^2*a*h*q - 1//2*T*a*b*h^3*q - 1//2*N*c1*h^3*q - 1//2*M*KT*h^3*q + 1//2*a*h^3*q)//(T*g*h), (-T^2*N*p - T*N*M*KN - T*N*f + T*N*g - T*N*h*p + T*alpha1 - N*M*KN*h - N*f*h + alpha1*h)//(T + h), (T^2*N*KN*c1*p + T*N*M*KN^2*c1 + T*N*KN*c1*f - T*N*KN*c1*g + T*N*KN*c1*h*p - T*KN*alpha1*c1 + N*M*KN^2*c1*h + N*KN*c1*f*h - KN*alpha1*c1*h)//(T + h), (-T^2*N*p - T*N*M*KN - T*N*f + T*N*g - T*N*h*p + T*alpha1 - N*M*KN*h - N*f*h + alpha1*h)//(T + h), (-T^2*N*c1*p - T*N*M*KN*c1 - T*N*c1*f + T*N*c1*g - T*N*c1*h*p + T*alpha1*c1 - N*M*KN*c1*h- N*c1*f*h + alpha1*c1*h)//(T + h), (T^2*N*a*b*p + T*N*M*KN*a*b + T*N*a*b*f - T*N*a*b*g + T*N*a*b*h*p - T*a*alpha1*b + N*M*KN*a*b*h + N*a*b*f*h - a*alpha1*b*h)//(T*N^2*q + N^2*h*q), (-T^2*N*p - T*N*M*KN - T*N*f + T*N*g - T*N*h*p + T*alpha1 - N*M*KN*h - N*f*h + alpha1*h)//(T + h), (-T^2*N*p - T*N*M*KN - T*N*f + T*N*g - T*N*h*p + T*alpha1 - N*M*KN*h - N*f*h + alpha1*h)//(T + h), (-T^2*N*p - T*N*M*KN - T*N*f + T*N*g - T*N*h*p + T*alpha1 - N*M*KN*h - N*f*h + alpha1*h)//(T + h), (T^2*N*KN*p + T*N*M*KN^2 + T*N*KN*f - T*N*KN*g + T*N*KN*h*p - T*KN*alpha1 + N*M*KN^2*h + N*KN*f*h - KN*alpha1*h)//(T + h), (-T^2*N*p - T*N*M*KN - T*N*f + T*N*g - T*N*h*p + T*alpha1 - N*M*KN*h - N*f*h + alpha1*h)//(T + h), (-T*L*I*q - T*L*gI*q + T*C*I*r2 + T*C*gI*r2 - L^2*I*ucte - L^2*gI*ucte - L*I*M*KL - L*I*m + L*I*pI- L*M*KL*gI - L*gI*m)//(I + gI), (-T*I*q - T*gI*q - L*I*KL - 2*L*I*ucte - L*KL*gI - 2*L*gI*ucte - I*M*KL - I*m + I*pI - M*KL*gI - gI*m)//(I + gI), (-T*L*I*q - T*L*gI*q + T*C*I*r2 + T*C*gI*r2 - L^2*I*ucte - L^2*gI*ucte - L*I*M*KL - L*I*m + L*I*pI - L*M*KL*gI - L*gI*m)//(I + gI), (-T*I*q - T*gI*q - L*I*KL - 2*L*I*ucte - L*KL*gI - 2*L*gI*ucte - I*M*KL - I*m + I*pI - M*KL*gI - gI*m)//(I + gI), (-T*L*I*q - T*L*gI*q + T*C*I*r2 + T*C*gI*r2 - L^2*I*ucte - L^2*gI*ucte - L*I*M*KL - L*I*m + L*I*pI - L*M*KL*gI - L*gI*m)//(I + gI), (-T*L*I*q*w - T*L*gI*q*w + T*C*I*r2*w + T*C*gI*r2*w - L^2*I*ucte*w - L^2*gI*ucte*w - L*I*M*KL*w - L*I*m*w + L*I*pI*w - L*M*KL*gI*w - L*gI*m*w)//(I + gI), (-T*L*I*q*w - T*L*gI*q*w + T*C*I*r2*w + T*C*gI*r2*w - L^2*I*ucte*w - L^2*gI*ucte*w - L*I*M*KL*w - L*I*m*w + L*I*pI*w - L*M*KL*gI*w - L*gI*m*w)//(I + gI), (-T*I*q - T*gI*q - L*I*KL - 2*L*I*ucte - L*KL*gI - 2*L*gI*ucte - I*M*KL - I*m + I*pI - M*KL*gI - gI*m)//(I + gI), (-T*I*q*w - T*gI*q*w - L*I*KL*w - 2*L*I*ucte*w - L*KL*gI*w - 2*L*gI*ucte*w - I*M*KL*w - I*m*w + I*pI*w - M*KL*gI*w - gI*m*w)//(I + gI), (-T*L*I*q - T*L*gI*q + T*C*I*r2 + T*C*gI*r2 - L^2*I*ucte - L^2*gI*ucte - L*I*M*KL - L*I*m + L*I*pI - L*M*KL*gI - L*gI*m)//(I + gI), (-T*I*q - T*gI*q - L*I*KL - 2*L*I*ucte - L*KL*gI - 2*L*gI*ucte - I*M*KL - I*m + I*pI - M*KL*gI - gI*m)//(I + gI), (-T*L*I*q - T*L*gI*q + T*C*I*r2 + T*C*gI*r2 - L^2*I*ucte - L^2*gI*ucte- L*I*M*KL - L*I*m + L*I*pI - L*M*KL*gI - L*gI*m)//(I + gI), (-T*I*q - T*gI*q - L*I*KL - 2*L*I*ucte - L*KL*gI - 2*L*gI*ucte - I*M*KL - I*m + I*pI - M*KL*gI - gI*m)//(I + gI), (-T*L*I*q - T*L*gI*q + T*C*I*r2 + T*C*gI*r2 - L^2*I*ucte - L^2*gI*ucte - L*I*M*KL - L*I*m + L*I*pI - L*M*KL*gI - L*gI*m)//(I + gI), (-T*I*q - T*gI*q - L*I*KL - 2*L*I*ucte - L*KL*gI - 2*L*gI*ucte - I*M*KL - I*m + I*pI - M*KL*gI - gI*m)//(I + gI)],
    )
)

# Examples arising from the wave-approach for MAPK

R, (KS00, KS01, KS10, FS01, FS10, FS11, K, F, S00, S01, S10, S11, a00, a01, a10, alpha01, alpha10, alpha11, b00, b01, b10, beta01, beta10, beta11, c0001, c0010, c0011, c0111, c1011, gamma0100, gamma1000, gamma1100, gamma1101, gamma1110) = polynomial_ring(QQ, ["KS00", "KS01", "KS10", "FS01", "FS10", "FS11", "K", "F", "S00", "S01", "S10", "S11", "a00", "a01", "a10", "alpha01", "alpha10", "alpha11", "b00", "b01", "b10", "beta01", "beta10", "beta11", "c0001", "c0010", "c0011", "c0111", "c1011", "gamma0100", "gamma1000", "gamma1100", "gamma1101", "gamma1110"])

push!(
    cases_simplification,
    Dict(
        :description => "New tide for MAPK, step 1",
        :gens => [S11, S00, F, S01 + S10, -FS10*beta10 - FS10*gamma1000 + F*S10*alpha10, -FS01*beta01 - FS01*gamma0100 + F*S01*alpha01],
    )
)

push!(
    cases_simplification,
    Dict(
        :description => "New tide for MAPK, step 2",
        :gens => [S11, S00, F, S01 + S10, FS10*beta10 + FS10*gamma1000 - F*S10*alpha10, FS01*beta01 + FS01*gamma0100 - F*S01*alpha01, FS01*beta01 + FS01*gamma0100 + FS10*beta10 + FS10*gamma1000 + FS11*beta11 + FS11*gamma1100 + FS11*gamma1101 + FS11*gamma1110 - F*S01*alpha01 - F*S10*alpha10 - F*S11*alpha11, FS01*beta01 + FS01*gamma0100 + FS10*beta10 + FS10*gamma1000 + FS11*beta11 + FS11*gamma1100 + FS11*gamma1101 + FS11*gamma1110 - F*S01*alpha01 - F*S10*alpha10 - F*S11*alpha11, KS00*b00 + KS00*c0001 + KS00*c0010 + KS00*c0011 + KS01*b01 + KS01*c0111 + KS10*b10 + KS10*c1011 - FS11*beta11 - FS11*gamma1100 - FS11*gamma1101 - FS11*gamma1110 - K*S00*a00 - K*S01*a01 - K*S10*a10 + F*S11*alpha11, KS00*b00 + KS00*c0001 + KS00*c0010 + KS00*c0011 + KS01*b01 + KS01*c0111 + KS10*b10 + KS10*c1011 - FS11*beta11 - FS11*gamma1100 - FS11*gamma1101 - FS11*gamma1110 - K*S00*a00 - K*S01*a01 - K*S10*a10 + F*S11*alpha11],
    )
)

push!(
    cases_simplification,
    Dict(
        :description => "New tide for MAPK, step 3",
        :gens => [ S11, S00, F, S01 + S10, FS10*beta10 + FS10*gamma1000 - F*S10*alpha10, FS01*beta01 + FS01*gamma0100 - F*S01*alpha01, FS11*beta11 + FS11*gamma1100 + FS11*gamma1101 + FS11*gamma1110 - F*S11*alpha11, KS00*b00 + KS00*c0001 + KS00*c0010 + KS00*c0011 + KS01*b01 + KS01*c0111 + KS10*b10 + KS10*c1011 - K*S00*a00 - K*S01*a01 - K*S10*a10, KS00*beta01*c0011 + KS00*c0011*gamma0100 + KS01*beta01*c0111 + KS01*c0111*gamma0100 + KS10*beta01*c1011 + KS10*c1011*gamma0100 + FS01*S01*alpha01*beta01 + FS01*S01*alpha01*gamma0100 + FS11*F*alpha01*beta11 + FS11*F*alpha01*gamma1100 + FS11*F*alpha01*gamma1101 + FS11*F*alpha01*gamma1110 + FS11*beta01*beta11 + FS11*beta11*gamma0100 - F^2*S11*alpha01*alpha11 - F*S01^2*alpha01^2 - F*S11*alpha11*beta01 - F*S11*alpha11*gamma0100, -KS00*beta01*c0011 - KS00*c0011*gamma0100 - KS01*beta01*c0111 - KS01*c0111*gamma0100 - KS10*beta01*c1011 - KS10*c1011*gamma0100 - FS01*S01*alpha01*beta01 - FS01*S01*alpha01*gamma0100 - FS11*F*alpha01*beta11 - FS11*F*alpha01*gamma1100 - FS11*F*alpha01*gamma1101 - FS11*F*alpha01*gamma1110 - FS11*beta01*beta11 - FS11*beta11*gamma0100 + F^2*S11*alpha01*alpha11 + F*S01^2*alpha01^2 + F*S11*alpha11*beta01 + F*S11*alpha11*gamma0100, KS00*beta01*c0011 + KS00*c0011*gamma0100 + KS01*beta01*c0111 + KS01*c0111*gamma0100 + KS10*beta01*c1011 + KS10*c1011*gamma0100 + FS01*S01*alpha01*beta01 + FS01*S01*alpha01*gamma0100 + FS11*F*alpha01*beta11 + FS11*F*alpha01*gamma1100 + FS11*F*alpha01*gamma1101 + FS11*F*alpha01*gamma1110 + FS11*beta01*beta11 + FS11*beta11*gamma0100 - F^2*S11*alpha01*alpha11 - F*S01^2*alpha01^2 - F*S11*alpha11*beta01 - F*S11*alpha11*gamma0100, -KS00*beta01*c0011 - KS00*c0011*gamma0100 - KS01*beta01*c0111 - KS01*c0111*gamma0100 - KS10*beta01*c1011 - KS10*c1011*gamma0100 - FS01*S01*alpha01*beta01 - FS01*S01*alpha01*gamma0100 - FS11*F*alpha01*beta11 - FS11*F*alpha01*gamma1100 - FS11*F*alpha01*gamma1101 - FS11*F*alpha01*gamma1110 - FS11*beta01*beta11 - FS11*beta11*gamma0100 + F^2*S11*alpha01*alpha11 + F*S01^2*alpha01^2 + F*S11*alpha11*beta01 + F*S11*alpha11*gamma0100, -KS00*beta01*c0011 - KS00*c0011*gamma0100 - KS01*beta01*c0111 - KS01*c0111*gamma0100 - KS10*beta01*c1011 - KS10*c1011*gamma0100 - FS01*S01*alpha01*beta01 - FS01*S01*alpha01*gamma0100 - FS11*F*alpha01*beta11 - FS11*F*alpha01*gamma1100 - FS11*F*alpha01*gamma1101 - FS11*F*alpha01*gamma1110 - FS11*beta01*beta11 - FS11*beta11*gamma0100 + F^2*S11*alpha01*alpha11 + F*S01^2*alpha01^2 + F*S11*alpha11*beta01 + F*S11*alpha11*gamma0100, -KS00*beta01*c0011 - KS00*c0011*gamma0100 - KS01*beta01*c0111 - KS01*c0111*gamma0100 - KS10*beta01*c1011 - KS10*c1011*gamma0100 - FS01*S01*alpha01*beta01 - FS01*S01*alpha01*gamma0100 - FS11*F*alpha01*beta11 - FS11*F*alpha01*gamma1100 - FS11*F*alpha01*gamma1101 - FS11*F*alpha01*gamma1110 - FS11*beta01*beta11 - FS11*beta11*gamma0100 + F^2*S11*alpha01*alpha11 + F*S01^2*alpha01^2 + F*S11*alpha11*beta01 + F*S11*alpha11*gamma0100 ],
    )
)

# Examples arising from the wave-approach for NFkB

R, (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, a1, a2, a3, c1, c1a, c1c, c2, c2a, c2c, c3, c3a, c3c, c4, c4a, c5, c5a, c6a, e1a, e2a, i1, i1a, k1, k2, k3, k_deg, k_prod, kv, t1, t2) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "a1", "a2", "a3", "c1", "c1a", "c1c", "c2", "c2a", "c2c", "c3", "c3a", "c3c", "c4", "c4a", "c5", "c5a", "c6a", "e1a", "e2a", "i1", "i1a", "k1", "k2", "k3", "k_deg", "k_prod", "kv", "t1", "t2"])

push!(
    cases_simplification,
    Dict(
        :description => "New tide for NFkB, step 1",
        :gens => [x12, x9, x7, x2, x10 + x13, x1 + x3, x6*i1*kv - x7*x11*a1, x1*k1 - x2*x8*k2, -x1*k1 + x2*x8*k2, x7*c1a - x12*c3a + c2a, -x7*c1 - x9*c3 + c2, -x1*x8*k1*k2 - x1*k1^2 + x2*x8^2*k2^2, x1*x8*k1*k2 + x1*k1^2 - x2*x8^2*k2^2],
    )
)

push!(
    cases_simplification,
    Dict(
        :description => "New tide for NFkB, step 2",
        :gens => [ x12, x9, x7, x2, x10 + x13, x1 + x3, x6*i1*kv - x7*x11*a1, x1*k1 - x2*x8*k2, x7*c1a - x12*c3a + c2a, x7*c1 + x9*c3 - c2, (x1*x8*k1*k2 + x1*k1^2 - x2*x8^2*k2^2)//x2, -x1*k_deg + x2*k3 - x3*k_deg + k_prod, (-x1^2*x8*k1^2*k2 - x1^2*k1^3 + x1*x2*x8^2*k1*k2^2 - x1*x2*k1^3)//x2^2, x1*x8^2*k1*k2^2 + x1*x8*k1^2*k2 + x1*k1^3 - x2*x8^3*k2^3, -x1*x8^2*k1*k2^2 - x1*x8*k1^2*k2 - x1*k1^3 + x2*x8^3*k2^3, -x1*x10*a2*k1 - x1*x13*a3*k1 + x2*x8*x10*a2*k2 + x2*x8*x13*a3*k2, x1*k1*k3 + x1*k1*k_deg - x2*x8*k2*k3 - x2*x8*k2*k_deg, x1*x8^2*k1*k2^2 + x1*x8*k1^2*k2 + x1*k1^3 - x2*x8^3*k2^3, -x1*k_deg + x2*k3 - x3*k_deg + k_prod, -x1*x10*a2*k1 - x1*x13*a3*k1 + x2*x8*x10*a2*k2 + x2*x8*x13*a3*k2, x1*k1*k3 + x1*k1*k_deg - x2*x8*k2*k3 - x2*x8*k2*k_deg, x6*c1a*i1*kv - x7*x11*a1*c1a - x7*c1a*c3a + x12*c3a^2 - c2a*c3a, x6*c1*i1*kv - x7*x11*a1*c1 - x7*c1*c3 - x9*c3^2 + c2*c3, x6*c1a*i1*kv - x7*x11*a1*c1a - x7*c1a*c3a + x12*c3a^2 - c2a*c3a, -x6*c1*i1*kv + x7*x11*a1*c1 + x7*c1*c3 + x9*c3^2 - c2*c3, x6*c1a*i1*kv - x7*x11*a1*c1a - x7*c1a*c3a + x12*c3a^2 - c2a*c3a, -x6*c1*i1*kv + x7*x11*a1*c1 + x7*c1*c3 + x9*c3^2 - c2*c3, x6*c1a*i1*kv - x7*x11*a1*c1a - x7*c1a*c3a + x12*c3a^2 - c2a*c3a, -x6*c1*i1*kv + x7*x11*a1*c1 + x7*c1*c3 + x9*c3^2 - c2*c3 ],
    )
)

push!(
    cases_simplification,
    Dict(
        :description => "New tide for NFkB, step 3",
        :gens => [ x12, x9, x7, x2, x8*k1*k2, k3 + k_deg, x10 + x13, x1 + x3, x10*a2 + x13*a3, x8*k2 + k1, x1*k1 + x2*k1, x6*i1*kv - x7*x11*a1, x7*c1a - x12*c3a + c2a, x7*c1 + x9*c3 - c2, x1*k_deg + x2*k_deg + x3*k_deg - k_prod, (x6*c1a*i1*kv - x7*x11*a1*c1a - x7*c1a*c3a + x12*c3a^2 - c2a*c3a)//(x6*i1*kv - x7*x11*a1), (x6*c1*i1*kv - x7*x11*a1*c1 - x7*c1*c3 - x9*c3^2 + c2*c3)//(x6*i1*kv - x7*x11*a1), -x8*c5*k1*k2 + x9*c4*k1*k2, -x8*c5*k2 + x9*c4*k2, -x1*x10*a2^2*k1 - x1*x13*a3^2*k1 + x2*x8*x10*a2^2*k2 + x2*x8*x13*a3^2*k2, -x1*x10*a2*k1*k_deg - x1*x13*a3*k1*k_deg + x2*x8*x10*a2*k2*k_deg + x2*x8*x13*a3*k2*k_deg, -x7*c1*c4*k1*k2 + x8*c5^2*k1*k2 - x9*c3*c4*k1*k2 - x9*c4*c5*k1*k2 + c2*c4*k1*k2, -x7*c1*c4*k2 + x8*c5^2*k2 - x9*c3*c4*k2 - x9*c4*c5*k2 + c2*c4*k2, -x2*x10*a2 - x2*x13*a3 - x2*k3 - x2*k_deg + x4*t1 + x5*t2, -x2*x10*a2 - x2*x13*a3 - x2*k3 - x2*k_deg + x4*t1 + x5*t2, -x2*x10*a2 - x2*x13*a3 - x2*k3 - x2*k_deg + x4*t1 + x5*t2, -x1*x10*a2*k1^2 - x1*x13*a3*k1^2 - x1*k1^2*k3 + x2*x8*x10*a2*k1*k2 + x2*x8*x13*a3*k1*k2 + x2*x8*k1*k2*k3 + x2*x8*k1*k2*k_deg, -x2*x10*a2 - x2*x13*a3 - x10*c5a - x10*i1a + x11*e1a + x12*c4a - x13*c6a + x14*e2a, -x1*k1*k_deg - x2*x10*a2*k1 - x2*x13*a3*k1 - x2*k1*k3 - x2*k1*k_deg + x4*k1*t1 + x5*k1*t2 + k1*k_prod, -x1*k_deg^2 - x2*x10*a2*k_deg - x2*x13*a3*k_deg - x2*k_deg^2 - x3*k_deg^2 + x4*k_deg*t1 + x5*k_deg*t2 + k_deg*k_prod, x2*x8*x10*a2*k1*k2 + x2*x8*x13*a3*k1*k2 + x2*x8*c5*k1*k2 + x2*x8*k1*k2*k3 + x2*x8*k1*k2*k_deg - x2*x9*c4*k1*k2 - x4*x8*k1*k2*t1 - x5*x8*k1*k2*t2, -x2*x10*a2 - x2*x13*a3 - x10*c5a - x10*i1a + x11*e1a + x12*c4a - x13*c6a + x14*e2a, -x2*x10*a2 - x2*x13*a3 - x10*c5a - x10*i1a + x11*e1a + x12*c4a - x13*c6a + x14*e2a, x5*i1*kv*t2 - x6*x10*a1*i1*kv - x6*x11*a1*i1*kv - x6*i1^2*kv + x7^2*x11*a1^2 - x7*x10*a1*i1a*kv + x7*x11^2*a1^2 + x7*x11*a1*e1a*kv + x13*c6a*i1*kv, x1*k_deg^2 - x2*x10*a2*k3 - x2*x13*a3*k3 - x2*k3^2 - 2*x2*k3*k_deg + x3*k_deg^2 + x4*k3*t1 + x5*k3*t2 - k_deg*k_prod, x1*k_deg^2 - x2*x10*a2*k3 - x2*x13*a3*k3 - x2*k3^2 - 2*x2*k3*k_deg + x3*k_deg^2 + x4*k3*t1 + x5*k3*t2 - k_deg*k_prod, x5*i1*kv*t2 - x6*x10*a1*i1*kv - x6*x11*a1*i1*kv - x6*i1^2*kv + x7^2*x11*a1^2 - x7*x10*a1*i1a*kv + x7*x11^2*a1^2 + x7*x11*a1*e1a*kv + x13*c6a*i1*kv, x1*k_deg^2 - x2*x10*a2*k3 - x2*x13*a3*k3 - x2*k3^2 - 2*x2*k3*k_deg + x3*k_deg^2 + x4*k3*t1 + x5*k3*t2 - k_deg*k_prod, x5*i1*kv*t2 - x6*x10*a1*i1*kv - x6*x11*a1*i1*kv - x6*i1^2*kv + x7^2*x11*a1^2 - x7*x10*a1*i1a*kv + x7*x11^2*a1^2 + x7*x11*a1*e1a*kv + x13*c6a*i1*kv, x5*i1*kv*t2 - x6*x10*a1*i1*kv - x6*x11*a1*i1*kv - x6*i1^2*kv + x7^2*x11*a1^2 - x7*x10*a1*i1a*kv + x7*x11^2*a1^2 + x7*x11*a1*e1a*kv + x13*c6a*i1*kv, x5*i1*kv*t2 - x6*x10*a1*i1*kv - x6*x11*a1*i1*kv - x6*i1^2*kv + x7^2*x11*a1^2 - x7*x10*a1*i1a*kv + x7*x11^2*a1^2 + x7*x11*a1*e1a*kv + x13*c6a*i1*kv, x1*k_deg^2 - x2*x10*a2*k3 - x2*x13*a3*k3 - x2*k3^2 - 2*x2*k3*k_deg + x3*k_deg^2 + x4*k3*t1 + x5*k3*t2 - k_deg*k_prod, x5*i1*kv*t2 - x6*x10*a1*i1*kv - x6*x11*a1*i1*kv - x6*i1^2*kv + x7^2*x11*a1^2 - x7*x10*a1*i1a*kv + x7*x11^2*a1^2 + x7*x11*a1*e1a*kv + x13*c6a*i1*kv, x1*k_deg^2 - x2*x10*a2*k3 - x2*x13*a3*k3 - x2*k3^2 - 2*x2*k3*k_deg + x3*k_deg^2 + x4*k3*t1 + x5*k3*t2 - k_deg*k_prod ],
    )
)


