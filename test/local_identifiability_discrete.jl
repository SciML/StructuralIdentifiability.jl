@testset "Discrete local identifiability, internal function" begin
 
    cases = []
    
    @parameters α β
    @variables t S(t) I(t) R(t) y(t)
    D = Difference(t; dt=1.0)
    
    eqs = [D(S) ~ S - β*S*I,
           D(I) ~ I + β*S*I - α*I,
           D(R) ~ R + α*I]
    @named sir = DiscreteSystem(eqs)
    push!(cases,
        Dict(
            :dds => sir,
            :res => Dict(S => true, I => true, R => false, α => true, β => true),
            :y => [y ~ I],
            :known_ic => Array{}[],
            :to_check => Array{}[]
        )
    )
    
    @parameters θ
    @variables t x(t) y(t)
    D = Difference(t; dt=1.0)
    
    eqs = [D(x) ~ θ*x^3]
    
    @named eqs = DiscreteSystem(eqs)
    push!(cases,
        Dict(
            :dds => eqs,
            :res => Dict(x => true, θ => true),
            :y => [y ~ x],
            :known_ic => Array{}[],
            :to_check => Array{}[]
        )
    )
    
    @parameters θ β
    @variables t x1(t) x2(t) y(t)
    D = Difference(t; dt=1.0)
    
    eqs = [D(x1) ~ x1 + x2,
           D(x2) ~ θ + β]
    
    @named eqs = DiscreteSystem(eqs)
    push!(cases,
        Dict(
            :dds => eqs,
            :res => Dict(x1 => true, x2 => true, θ => false, β => false),
            :y => [y ~ x1],
            :known_ic => Array{}[],
            :to_check => Array{}[]
        )
    )
    
    @parameters a b c d
    @variables t x1(t) x2(t) u(t)
    D = Difference(t; dt=1.0)
    
    eqs = [
        D(x1) ~ a * x1 - b * x1 * x2 + u,
        D(x2) ~ -c * x2 + d * x1 * x2
    ]
    
    @named lv = DiscreteSystem(eqs)
    push!(cases,
        Dict(
            :dds => lv,
            :res => Dict(a => true, b => false, c => true, d => true, x1 => true, x2 => false),
            :y => [y ~ x1],
            :known_ic => Array{}[],
            :to_check => Array{}[]
        )
    )
    
    push!(cases,
        Dict(
            :dds => lv,
            :res => Dict(b * x2 => true),
            :y => [y ~ x1],
            :known_ic => Array{}[],
            :to_check => [b * x2]
        )
    )
    
    push!(cases,
        Dict(
            :dds => lv,
            :res => Dict(a => true, b => true, c => true, d => true, x1 => true, x2 => true),
            :y => [y ~ x1],
            :known_ic => [x2],
            :to_check => Array{}[]
        )
    )
    
    # Example 1 from https://doi.org/10.1016/j.automatica.2008.03.019
    @parameters theta1 theta2
    @variables t x1(t) x2(t) u(t) y(t)
    D = Difference(t; dt=1.0)
    
    eqs = [
        D(x1) ~ (1 + theta1) * x1 + x2,
        D(x2) ~ (1 - theta2) * x1 + x2^2 + u
    ]
    
    @named abmd1 = DiscreteSystem(eqs)
    push!(cases,
        Dict(
            :dds => abmd1,
            :res => Dict(x1 => true, x2 => true, theta1 => true, theta2 => true),
            :y => [y ~ x1],
            :known_ic => Array{}[],
            :to_check => Array{}[]
        )
    )
    
    # Example 2 from https://doi.org/10.1016/j.automatica.2008.03.019
    @parameters theta1 theta2 theta3
    @variables t x1(t) x2(t) u(t) y(t)
    D = Difference(t; dt=1.0)
    
    eqs = [
        D(x1) ~ theta1 * x1^2 + theta2 * x2 + u,
        D(x2) ~ theta3 * x1
    ]
    
    @named abmd2 = DiscreteSystem(eqs)
    push!(cases,
        Dict(
            :dds => abmd2,
            :res => Dict(theta1 => true, theta2 => false, theta3 => false, x1 => true, x2 => false),
            :y => [y ~ x1],
            :known_ic => Array{}[],
            :to_check => Array{}[]
        )
    )
    
    for c in cases
        @test assess_local_identifiability(c[:dds]; measured_quantities=c[:y], known_ic=c[:known_ic], funcs_to_check=c[:to_check]) == c[:res]
    end
end
