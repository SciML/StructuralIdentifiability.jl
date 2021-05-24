read "../IdentifiabilityODE.mpl";

sys := [
diff(EGF_EGFR(t), t) = reaction_1_k1*EGF_EGFR(t) - reaction_9_k1*EGF_EGFR(t) - reaction_1_k2*EGF_EGFR(t),
diff(EGFR(t), t) = -reaction_1_k1*EGF_EGFR(t) + reaction_1_k2*EGF_EGFR(t) - EGFR(t)*EGFR_turnover + EGFR_turnover*pro_EGFR(t),
diff(pEGFR(t), t) = -reaction_4_k1*pEGFR(t) + reaction_9_k1*EGF_EGFR(t) - pEGFR(t)*Akt(t)*reaction_2_k1 + reaction_3_k1*pEGFR_Akt(t) + pEGFR_Akt(t)*reaction_2_k2,
diff(pAkt_S6(t), t) = pAkt(t)*reaction_5_k1*S6(t) - reaction_6_k1*pAkt_S6(t) - pAkt_S6(t)*reaction_5_k2,
diff(pEGFR_Akt(t), t) = pEGFR(t)*Akt(t)*reaction_2_k1 - reaction_3_k1*pEGFR_Akt(t) - pEGFR_Akt(t)*reaction_2_k2,
diff(pAkt(t), t) = -pAkt(t)*reaction_7_k1 - pAkt(t)*reaction_5_k1*S6(t) + reaction_6_k1*pAkt_S6(t) + reaction_3_k1*pEGFR_Akt(t) + pAkt_S6(t)*reaction_5_k2,
diff(S6(t), t) = pS6(t)*reaction_8_k1 - pAkt(t)*reaction_5_k1*S6(t) + pAkt_S6(t)*reaction_5_k2,
diff(pS6(t), t) = -pS6(t)*reaction_8_k1 + reaction_6_k1*pAkt_S6(t),
diff(Akt(t), t) = pAkt(t)*reaction_7_k1 - pEGFR(t)*Akt(t)*reaction_2_k1 + pEGFR_Akt(t)*reaction_2_k2,
y1(t) = pEGFR(t)*a1 + a1*pEGFR_Akt(t),
y2(t) = a2*pAkt(t) + a2*pAkt_S6(t),
y3(t) = pS6(t)*a3
];
CodeTools[CPUInfo](IdentifiabilityODE(sys, GetParameters(sys)));