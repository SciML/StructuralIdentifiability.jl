using Dates
using Logging
using Oscar

struct ODE
    poly_ring
    x_vars
    u_vars
    parameters
    equations
    
    function ODE(eqs, inputs, field = QQ)
        #Initialize ODE
        #equations is a dictionary x_i => f_i(x, u, params)

        poly_ring = parent(collect(values(eqs))[1])
        x_vars = collect(keys(eqs))
        u_vars = inputs
        parameters = filter(v -> (!(v in x_vars) && !(v in u_vars)), gens(poly_ring))
        new(poly_ring, x_vars, u_vars, parameters, eqs)
    end
end                
