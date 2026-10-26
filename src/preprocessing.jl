function preprocess_ode(ode::ODE{P}; output_saturation = true) where {P <: MPolyRingElem}
    submodels = find_submodels(ode)
    if length(submodels) > 0
        @info "Note: the input model has nontrivial submodels. If the computation for the full model will be too heavy, you may want to try to first analyze one of the submodels. They can be produced using function `find_submodels`"
    end

    new_ode = ode

    if length(ode.u_vars) > 0 && output_saturation
        @info "Model contains outputs, we preprocess it with output saturation. This typically speeds up computations, in the case you want to turn this off, change `output_saturation` option"
        saturation_orders = propose_orders(ode)
        new_ode = saturate_outputs(ode, saturation_orders)
        @debug "Added $(length(new_ode.y_vars) - length(ode.y_vars)) new outputs"
    end

    return new_ode
end
