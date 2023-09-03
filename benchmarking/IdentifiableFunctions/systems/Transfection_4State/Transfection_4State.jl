# Transfection_4State
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	mRNAenz'(t) = enz(t)*d2*mRNA(t) - mRNAenz(t)*d3,
	enz'(t) = -enz(t)*d2*mRNA(t) + mRNAenz(t)*d3,
	GFP'(t) = -b*GFP(t) + kTL*mRNA(t),
	mRNA'(t) = -enz(t)*d2*mRNA(t) - mRNA(t)*d1,
	y1(t) = GFP(t)
)
