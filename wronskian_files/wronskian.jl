using Dates
using IterTools
using Logging
using Oscar

include("../util.jl")
include("../elimination.jl")
include("../ODE.jl")
include("../io_equation.jl")

#=
We want to construct the Wronskian, translate the Mathematica code into Julia
=#

function monomial_compress(ioequation::fmpq_mpoly,ode::ODE)

	params = ode.parameters # Extract the list of parameters from the io equation
	variables = gens(parent(ioequation))[findall(!in(map(string,params)),map(string,gens(parent(ioequation))))]
	coeffdict = extract_coefficients(ioequation,variables)
	expvect = collect(keys(coeffdict))
	coeffs = collect(values(coeffdict))
	terms = map(x->prod(variables.^x),expvect)
	# parameterRing, parametersymbols = PolynomialRing(base_ring(parent(ioequation)),map((x)->string(x),params)) # Create a ring of parameters
	# refactored_ring, variable_symbols = PolynomialRing(parameterRing,map((x)->string(x),vcat(ode.x_vars,ode.u_vars)))
	# print(refactored_ring) # the refactored ring is correct, but it breaks parent_ring_change
	# refactored_ioequation = parent_ring_change(ioequation,refactored_ring ) # change the input io equation so that it lives in the ring of variables over the ring of parameters over the base ring

	# We want to factor each monomial in the form f(p)*g(x,u) where overall numerical factors are in g.
	# terms = collect(monomials(refactored_ioequation)) # These are the g(x,u) above
	# coeffs = collect(coefficients(refectored_ioequation)) # These are the f(p) above (At the moment the common factors are here rather than in the terms list.
	consts = map((x)->collect(coefficients(x)),coeffs) # Take the list of constants from the coefficients
	gcds = map((x)->reduce(gcd,x),consts) # take the gcds of each of these lists to get the overall numerical factors
	reducedcoeffs = coeffs.//gcds # divide the coefficients by the overall numerical factors
	modifiedterms = terms.*gcds # multiply the terms by the overall numerical factors
	
	# We have now gotten each monomial into the form f(p)*g(x,u).


	#=
	as a sanity check, taking the inner product of the terms list and the coeffs list should yield the original io equation
	=#


	# Now, we want to loop through the coefficients, and if we come across a repeat, we want to add the corresponding terms and delete the repeat.
	foundDuplicateIndices = Vector{Int}()
	collectedterms = Vector()
	condensed_coefficients = Vector()
	for i = 1:length(reducedcoeffs)
		if !(i in foundDuplicateIndices)
			testcoeff = reducedcoeffs[i]
			locations = findall(x->x==testcoeff,reducedcoeffs)
			append!(foundDuplicateIndices,locations)
			condensed_term = sum(modifiedterms[locations])
			push!(collectedterms,condensed_term)
			push!(condensed_coefficients,testcoeff)

		end
	end
	return (condensed_coefficients, collectedterms)
end

#----------------------------------------------------------------------------------------------------

function wronskian()

	return

end
