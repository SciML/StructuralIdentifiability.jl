# Bruno2016
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	OHbetaio'(t) = cry(t)*kcryOH + zea(t)*kzea + OHbeta10(t)*kOHbeta10,
	betaio'(t) = cry(t)*kcrybeta + beta10(t)*kbeta10 + kbeta*beta(t),
	cry'(t) = -cry(t)*kcrybeta - cry(t)*kcryOH,
	OHbeta10'(t) = cry(t)*kcrybeta + zea(t)*kzea - OHbeta10(t)*kOHbeta10,
	beta10'(t) = cry(t)*kcryOH - beta10(t)*kbeta10 + kbeta*beta(t),
	zea'(t) = -zea(t)*kzea,
	beta'(t) = -kbeta*beta(t),
	y1(t) = beta(t),
	y2(t) = beta10(t)
)
