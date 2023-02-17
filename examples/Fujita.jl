# K. A. Fujita, Y. Toyoshima, S. Uda, Y. I. Ozaki, H. Kubota, and S. Kuroda
# Decoupling of receptor and downstream signals in the Akt pathway by its low-pass filter characteristics
# https://doi.org/10.1126/scisignal.2000810

using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

#! format: off
ode = @ODEmodel(
    EGFR'(t) = EGFR_turnover * pro_EGFR(t) + EGF_EGFR(t) * reaction_1_k2 - EGFR(t) * EGFR_turnover - EGF_EGFR(t) * reaction_1_k1,
    pEGFR'(t) = EGF_EGFR(t) * reaction_9_k1 - pEGFR(t) * reaction_4_k1 + pEGFR_Akt(t) * reaction_2_k2 + pEGFR_Akt(t) * reaction_3_k1 - Akt(t) * pEGFR(t) * reaction_2_k1,
    pEGFR_Akt'(t) = Akt(t) * pEGFR(t) * reaction_2_k1 - pEGFR_Akt(t) * reaction_3_k1 - pEGFR_Akt(t) * reaction_2_k2,
    Akt'(t) = pAkt(t) * reaction_7_k1 + pEGFR_Akt(t) * reaction_2_k2 - Akt(t) * pEGFR(t) * reaction_2_k1,
    pAkt'(t) = pAkt_S6(t) * reaction_5_k2 - pAkt(t) * reaction_7_k1 + pAkt_S6(t) * reaction_6_k1 + pEGFR_Akt(t) * reaction_3_k1 - S6(t) * pAkt(t) * reaction_5_k1,
    S6'(t) = pAkt_S6(t) * reaction_5_k2 + pS6(t) * reaction_8_k1 - S6(t) * pAkt(t) * reaction_5_k1,
    pAkt_S6'(t) = S6(t) * pAkt(t) * reaction_5_k1 - pAkt_S6(t) * reaction_6_k1 - pAkt_S6(t) * reaction_5_k2,
    pS6'(t) = pAkt_S6(t) * reaction_6_k1 - pS6(t) * reaction_8_k1,
    EGF_EGFR'(t) = EGF_EGFR(t) * reaction_1_k1 - EGF_EGFR(t) * reaction_9_k1 - EGF_EGFR(t) * reaction_1_k2,
    y1(t) = a1 * (pEGFR(t) + pEGFR_Akt(t)),
    y2(t) = a2 * (pAkt(t) + pAkt_S6(t)),
    y3(t) = a3 * pS6(t),
)
#! format: on

@time println(assess_identifiability(ode))
