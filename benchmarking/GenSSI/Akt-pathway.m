function model = SMTH()
syms EGF_EGFR EGFR pEGFR pAkt_S6 pEGFR_Akt pAkt S6 pS6 Akt
syms a2 a3 reaction_1_k1 reaction_4_k1 reaction_9_k1 reaction_1_k2 reaction_7_k1 reaction_6_k1 a1 reaction_3_k1 reaction_5_k1 reaction_2_k1 EGFR_turnover reaction_2_k2 reaction_5_k2 reaction_8_k1
syms EGF_EGFR0 EGFR0 pEGFR0 pAkt_S60 pEGFR_Akt0 pAkt0 S60 pS60 Akt0
syms pro_EGFR
model.sym.p = [a2; a3; reaction_1_k1; reaction_4_k1; reaction_9_k1; reaction_1_k2; reaction_7_k1; reaction_6_k1; a1; reaction_3_k1; reaction_5_k1; reaction_2_k1; EGFR_turnover; reaction_2_k2; reaction_5_k2; reaction_8_k1; EGF_EGFR0; EGFR0; pEGFR0; pAkt_S60; pEGFR_Akt0; pAkt0; S60; pS60; Akt0];
model.sym.x = [EGF_EGFR; EGFR; pEGFR; pAkt_S6; pEGFR_Akt; pAkt; S6; pS6; Akt];
model.sym.g = [pro_EGFR];
model.sym.x0 = [EGF_EGFR0; EGFR0; pEGFR0; pAkt_S60; pEGFR_Akt0; pAkt0; S60; pS60; Akt0];
model.sym.xdot = [reaction_1_k1*EGF_EGFR - reaction_9_k1*EGF_EGFR - reaction_1_k2*EGF_EGFR
-reaction_1_k1*EGF_EGFR + reaction_1_k2*EGF_EGFR - EGFR*EGFR_turnover + EGFR_turnover*pro_EGFR
-reaction_4_k1*pEGFR + reaction_9_k1*EGF_EGFR - pEGFR*Akt*reaction_2_k1 + reaction_3_k1*pEGFR_Akt + pEGFR_Akt*reaction_2_k2
pAkt*reaction_5_k1*S6 - reaction_6_k1*pAkt_S6 - pAkt_S6*reaction_5_k2
pEGFR*Akt*reaction_2_k1 - reaction_3_k1*pEGFR_Akt - pEGFR_Akt*reaction_2_k2
-pAkt*reaction_7_k1 - pAkt*reaction_5_k1*S6 + reaction_6_k1*pAkt_S6 + reaction_3_k1*pEGFR_Akt + pAkt_S6*reaction_5_k2
pS6*reaction_8_k1 - pAkt*reaction_5_k1*S6 + pAkt_S6*reaction_5_k2
-pS6*reaction_8_k1 + reaction_6_k1*pAkt_S6
pAkt*reaction_7_k1 - pEGFR*Akt*reaction_2_k1 + pEGFR_Akt*reaction_2_k2];
model.sym.y = [pEGFR*a1 + a1*pEGFR_Akt
a2*pAkt + a2*pAkt_S6
pS6*a3];
end