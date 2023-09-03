# LeukaemiaLeon2021
#! format: off
using StructuralIdentifiability

system = @ODEmodel(
	P'(t) = (2*ap*rhop*taop*P(t) - ks*rhop*taop*P(t)^2 - ks*rhop*taop*P(t)*I(t) - ks*P(t)^2 - ks*P(t)*I(t) - rhop*taop*P(t) - P(t))//(ks*taop*P(t) + ks*taop*I(t) + taop),
	C'(t) = (B(t)*rhoc*taoc*C(t) + rhoc*taoc*L(t)*C(t) + I(t)*taoc*rhob*C(t) - C(t))//taoc,
	L'(t) = -alpha*L(t)*C(t) + rhol*L(t),
	I'(t) = (-alpha*taoi*ks*taop*P(t)*I(t)*beta*C(t) - alpha*taoi*ks*taop*I(t)^2*beta*C(t) - alpha*taoi*taop*I(t)*beta*C(t) - taoi*ks*taop*P(t)*I(t)*rhoi - taoi*ks*taop*I(t)^2*rhoi + taoi*ks*P(t)^2 + taoi*ks*P(t)*I(t) + 2*taoi*ai*taop*I(t)*rhoi - taoi*taop*I(t)*rhoi + taoi*P(t) - ks*taop*P(t)*I(t) - ks*taop*I(t)^2 - taop*I(t))//(taoi*ks*taop*P(t) + taoi*ks*taop*I(t) + taoi*taop),
	B'(t) = (-alpha*B(t)*taoi*taob*C(t) - B(t)*taoi + taob*I(t))//(taoi*taob),
	y1(t) = C(t),
	y2(t) = L(t),
	y3(t) = B(t) + P(t) + I(t)
)
