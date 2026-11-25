const PROJECTION_VARNAME = "rand_proj_var"

"""
    find_extra_projection(ode, io_projections)

Finds an extra IO-projection to restore the primality of the ideal, see Remark 3 in https://arxiv.org/pdf/2111.00991.
The code is based on Algorithm 4.5 in the paper.
Input:
- `ode` - the ODE system
- `io_projections` - a projection-based representation of the ODE as returned by `find_ioequations`

Output: a dictionary extending `io_projections` by an extra projection
"""

function find_extra_projection(ode::ODE{P}, io_projections::Dict{P, P}) where {P <: MPolyRingElem{<:FieldElem}}
    ring = parent(first(values(io_projections)))
    new_ring, _ = polynomial_ring(Nemo.QQ, [map(var_to_str, gens(ring))..., PROJECTION_VARNAME])
    proj_var = gens(new_ring)[end]
    new_io_projections = Dict(
        parent_ring_change(k, new_ring) => parent_ring_change(v, new_ring)
            for (k, v) in io_projections
    )
    gpg = generic_point_generator = ODEPointGenerator{P}(ode, new_ring)

    sampling_range = 5
    while true
        nonlinear_leaders = [v for (v, eq) in new_io_projections if degree(eq, v) > 1]
        @assert !isempty(nonlinear_leaders)
        extra_projection = sum(rand(1:sampling_range) * v for v in nonlinear_leaders)
        @debug "Extra projections: $extra_projection"
        gpg = generator_var_change(gpg, proj_var, extra_projection, one(new_ring))

        projection_equation = proj_var - extra_projection

        for v in nonlinear_leaders
            projection_equation = eliminate_var(projection_equation, new_io_projections[v], v, gpg)
        end
        @debug "Check primality"
        projection_pullback = evaluate(projection_equation, [proj_var], [extra_projection])
        if check_primality(new_io_projections, [projection_pullback])
            @debug "Single component of highest dimension isolated, returning"
            new_io_projections[proj_var] = projection_equation
            return new_io_projections
        end
        sampling_range = 2 * sampling_range
    end
    return
end
