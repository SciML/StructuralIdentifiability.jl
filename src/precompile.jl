
@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    using Logging
    using ModelingToolkit
    @parameters a01 a21 a12
    @variables t x0(t) x1(t) y(t)
    D = Differential(t)
    eqs = [D(x0) ~ -(a01 + a21) * x0 + a12 * x1, D(x1) ~ a21 * x0 - a12 * x1]
    de = ODESystem(eqs, t, name = :Test)
    @compile_workload begin
        with_logger(Logging.ConsoleLogger(Logging.Warn)) do
            # all calls in this block will be precompiled, regardless of whether
            # they belong to your package or not (on Julia 1.8 and higher)
            restart_logging(loglevel = Logging.Warn)
            ode = @ODEmodel(
                x1'(t) = -(a01 + a21) * x1(t) + a12 * x2(t) + u(t),
                x2'(t) = a21 * x1(t) - a12 * x2(t) - x3(t) / b,
                x3'(t) = x3(t),
                y(t) = x2(t)
            )
            assess_identifiability(ode)
            assess_identifiability(de; measured_quantities = [x0])
            assess_identifiability(de; measured_quantities = [y ~ x0])
            find_identifiable_functions(ode, with_states = true)
        end
    end
end
