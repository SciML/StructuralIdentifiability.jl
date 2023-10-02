# Transfection_4State
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	mRNA'(t) = -enz(t)*d2*mRNA(t) - mRNA(t)*d1,
	GFP'(t) = -b*GFP(t) + kTL*mRNA(t),
	enz'(t) = -enz(t)*d2*mRNA(t) + mRNAenz(t)*d3,
	mRNAenz'(t) = enz(t)*d2*mRNA(t) - mRNAenz(t)*d3,
	y1(t) = GFP(t)
)
