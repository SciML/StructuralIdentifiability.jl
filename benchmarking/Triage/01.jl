# https://doi.org/10.1051/meca/2015054
# eq. 3-5

using StructuralIdentifiability

ode = @ODEmodel(
    Ti'(t) = ((Tap(t) - Ti(t))/Rap - (Ti(t) - Te(t))/Rth)/Cth,

    Tap'(t) = (Pelec(t) - (Tap(t) - Ti(t))/Rap)/Cap,

    Te'(t) = 0,
    Pelec'(t) = 0,

    y(t) = Ti(t)
)

find_identifiable_functions(
    ode;
    known_ic = [Te, Pelec],
    with_states = true
)
