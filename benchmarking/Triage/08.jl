# https://doi.org/10.3390/en4040582

using StructuralIdentifiability

ode = @ODEmodel(
           Ud'(t) = a*u(t),
           Upn'(t) = -Upn(t)/(Rpn*Cpn) + u(t)/Cpn,
           y(t) = Uoc - Ud(t) - Upn(t) - Ro*u(t)
       )

id = find_identifiable_functions(ode; with_states=true)
