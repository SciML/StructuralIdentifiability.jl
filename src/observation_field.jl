"""
    observation_field(ode::ODE; options...)

Find simple generators of the observation field of the given ODE system.

## Options

This functions takes the following optional arguments:
- `known_ic`: a list of functions whose initial conditions are assumed to be known,
  then the returned identifiable functions will be functions of parameters and
  initial conditions, not states (this is an experimental functionality).
- `prob_threshold`: A float in the range from 0 to 1, the probability of correctness. Default
  is `0.99`.
- `seed`: The rng seed. Default value is `42`.
- `loglevel` - the minimal level of log messages to display (`Logging.Info` by default)

## Example

```jldoctest
using StructuralIdentifiability

ode = @ODEmodel(
    x0'(t) = -(a01 + a21) * x0(t) + a12 * x1(t),
    x1'(t) = a21 * x0(t) - a12 * x1(t),
    y(t) = x0(t)
);

observation_field(ode)

# Output

4-element Vector{AbstractAlgebra.Generic.FracFieldElem{Nemo.QQMPolyRingElem}}:
 x0(t)
 a01 + a12 + a21
 a01*a12
 x0(t)*a12 + x1(t)*a12
```

"""
function observation_field(
        ode::ODE{T};
        known_ic::Vector{<:ExtendedFraction{T}} = Vector{ExtendedFraction{T}}(),
        prob_threshold::Float64 = 0.99,
        seed = 42,
        loglevel = Logging.Info,
    ) where {T <: MPolyRingElem{QQFieldElem}}
    return find_identifiable_functions(ode, with_states = true, known_ic = known_ic, prob_threshold = prob_threshold, seed = seed, loglevel = loglevel)
end
