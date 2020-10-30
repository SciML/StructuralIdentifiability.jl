include("util.jl")

#------------------------------------------------------------------------------

function truncate_matrix(M::MatElem{<: Generic.AbsSeriesElem}, prec::Int)
    return map(e -> truncate(e, prec), M)
end

#------------------------------------------------------------------------------

function matrix_set_prec!(M::MatElem{<: Generic.AbsSeriesElem}, prec::Int)
    map(e -> set_prec!(e, prec), M)
end

#------------------------------------------------------------------------------

function ps_matrix_const_term(M::MatElem{<: Generic.AbsSeriesElem})
    return map(e -> coeff(e, 0), M)
end

#------------------------------------------------------------------------------

function _matrix_inv_newton_iteration(M::MatElem{T}, Minv::MatElem{T}) where T <: Generic.AbsSeriesElem{<: Generic.FieldElem}
    """
    Performs a single step of Newton iteration for inverting M with 
    Minv being a partial result
    """
    return 2 * Minv - Minv * M * Minv
end

#--------------------------------------

function ps_matrix_inv(M::MatElem{<: Generic.AbsSeriesElem{<: Generic.FieldElem}}, prec::Int=-1)
    """
    Input:
        - M - a square matrix with entries in a univariate power series ring
          it is assumed that M(0) is invertible and all entries having the same precision
    Output: the inverse of M computed up to prec (the precision of the matrix if nothing)
    """
    const_term = ps_matrix_const_term(M)
    prec = (prec == -1) ? precision(M[1, 1]) : prec
    power_series_ring = base_ring(parent(M))
    result = map(a -> power_series_ring(a), Base.inv(const_term))
    cur_prec = 1
    while cur_prec < prec
        result = _matrix_inv_newton_iteration(M, result)
        cur_prec *= 2
    end
    return result
end

#------------------------------------------------------------------------------

function ps_diff(ps::Generic.AbsSeriesElem{<: Generic.RingElem})
    """
    Input: ps - (absolute capped) unvariate power series
    Output: the derivative of ps
    """
    result = zero(parent(ps))
    set_prec!(result, precision(ps))
    for exp in 1:(precision(ps) - 1)
        setcoeff!(result, exp - 1, coeff(ps, exp) * exp)
    end
    return result
end

#------------------------------------------------------------------------------

function ps_integrate(ps::Generic.AbsSeriesElem{<: Generic.FieldElem})
    """
    Input: ps - (absolute capped) unvariate power series
    Output: the integral of ps without constant term
    """
    result = zero(parent(ps))
    set_prec!(result, precision(ps) + 1)
    for exp in 0:(precision(ps) - 1)
        setcoeff!(result, exp + 1, coeff(ps, exp) // (exp + 1))
    end
    return result
end

#------------------------------------------------------------------------------

function ps_matrix_log(M::MatElem{<: Generic.AbsSeriesElem{<: Generic.FieldElem}})
    """
    Input:
        - M - a square matrix with entries in a univariate power series ring
          it is assumed that M(0) is the identity
    Output: the natural log of M
    """
    const_term = ps_matrix_const_term(M)
    if const_term != one(parent(const_term))
       throw(Base.DomainError("Constant term must be the identity matrix"))
    end
    return map(
        e -> ps_integrate(e),
        map(e -> ps_diff(e), M - one(parent(M))) *
        ps_matrix_inv(M)
    )
end

#------------------------------------------------------------------------------

function _matrix_homlinear_de_newton_iteration(
        A::MatElem{<: Generic.AbsSeriesElem{T}},
        Y::MatElem{<: Generic.AbsSeriesElem{T}},
        Z::MatElem{<: Generic.AbsSeriesElem{T}},
        cur_prec::Int
    ) where T <: Generic.FieldElem
    Yprime = map(ps_diff, Y)
    Z = Z + truncate_matrix(Z * (one(parent(A)) - Y * Z), cur_prec)
    Y = Y - truncate_matrix(Y * map(ps_integrate, Z * (Yprime - truncate_matrix(A, 2 * cur_prec - 1) * Y)), 2 * cur_prec)
    return (Y, Z)
end

#--------------------------------------

function ps_matrix_homlinear_de(
        A::MatElem{<: Generic.AbsSeriesElem{T}},
        Y0::MatElem{<: T},
        prec::Int=-1
    ) where T <: Generic.FieldElem
    """
    Input:
        - A - a square matrix with entries in a univariate power series ring
        - Y0 - a square invertible matrix over the base field
    Output: matrix Y such that Y' = AY up to precision of A - 1 and Y(0) = Y0
    """
    prec = (prec == -1) ? precision(A[1, 1]) : prec
    ps_ring = base_ring(parent(A))
    cur_prec = 1
    Y = (one(parent(A)) + gen(ps_ring) * truncate_matrix(A, cur_prec)) * map(e -> truncate(ps_ring(e), cur_prec), Y0)
    Z = map(e -> truncate(ps_ring(e), cur_prec), Base.inv(Y0))
    while cur_prec < prec
        matrix_set_prec!(Y, 2 * cur_prec)
        matrix_set_prec!(Z, cur_prec)
        Y, Z = _matrix_homlinear_de_newton_iteration(A, Y, Z, cur_prec)
        cur_prec *= 2
    end
    return Y, Z
end

#------------------------------------------------------------------------------

function _variation_of_constants(
        A::MatElem{<: Generic.AbsSeriesElem{T}},
        B::MatElem{<: Generic.AbsSeriesElem{T}},
        Yh::MatElem{<: Generic.AbsSeriesElem{T}},
        Zh::MatElem{<: Generic.AbsSeriesElem{T}},
        Y0::MatElem{<: T},
        prec::Int
    ) where T <: Generic.FieldElem
    Zh += truncate_matrix(Zh * (one(parent(A)) - Yh * Zh), prec)
    Y_particular = truncate_matrix(Yh * map(ps_integrate, Zh * B), prec)
    return Y_particular + Yh * map(e -> parent(A[1, 1])(e), Y0)
end

#--------------------------------------

function ps_matrix_linear_de(
        A::MatElem{<: Generic.AbsSeriesElem{T}},
        B::MatElem{<: Generic.AbsSeriesElem{T}},
        Y0::MatElem{<: T},
        prec::Int=-1
    ) where T <: Generic.FieldElem
    """
    Input:
        - A, B - square matrices with entries in a univariate power series ring
        - Y0 - a matrix over the base field with the rows number the same as A
    Output: matrix Y such that Y' = AY + B up to precision of A - 1 and Y(0) = Y0
    """
    prec = (prec == -1) ? precision(A[1, 1]) : prec
    n = nrows(A)
    identity = one(MatrixSpace(base_ring(parent(Y0)), n, n))
    Yh, Zh = ps_matrix_homlinear_de(A, identity, prec)
    matrix_set_prec!(Zh, prec)
    return _variation_of_constants(A, B, Yh, Zh, Y0, prec)
end

#------------------------------------------------------------------------------

function ps_ode_solution(
        equations::Array{P, 1},
        ic::Dict{P, T},
        inputs::Dict{P, Array{T, 1}},
        prec::Int
    ) where {T <: Generic.FieldElem, P <: MPolyElem{T}}
    """
    Input
        - equations - a system of the form A(x, u, mu)x' - B(x, u, mu) = 0,
                      where A is a generically nonsingular square matrix
        - ic - initial conditions for x's (dictionary)
        - inputs - power series for inputs represented as arrays (dictionary)
        - prec - precision of the solution
        Assumption: A is nonzero at zero
    Output: power series solution of the system
    """
    n = length(equations)
    ring = parent(equations[1])
    S = MatrixSpace(ring, n, n)
    Sv = MatrixSpace(ring, n, 1)
    Svconst = MatrixSpace(base_ring(ring), n, 1)
    eqs = Sv(equations)
    
    x_vars = filter(v -> ("$(v)_dot" in map(string, gens(ring))), gens(ring))
    x_vars = [x for x in x_vars]
    x_dot_vars = [str_to_var(var_to_str(x) * "_dot", ring) for x in x_vars]

    Jac_dots = S([derivative(p, xd) for p in equations, xd in x_dot_vars])
    Jac_xs = S([derivative(p, x) for p in equations, x in x_vars])

    ps_ring, t = PowerSeriesRing(base_ring(ring), prec, "t"; model=:capped_absolute)
    solution = Dict()
    for (u, coeffs) in inputs
        solution[u] = sum([coeffs[i] * t^(i - 1) for i in 1:length(coeffs)])
    end
    for x in x_vars
        solution[x] = ps_ring(ic[x])
    end
    for xd in x_dot_vars
        solution[xd] = ps_ring(0)
        set_prec!(solution[xd], 1)
    end

    cur_prec = 1
    while cur_prec < prec
        for i in 1:length(x_vars)
            set_prec!(solution[x_vars[i]], 2 * cur_prec)
            set_prec!(solution[x_dot_vars[i]], 2 * cur_prec)
        end
        eval_point = [solution[v] for v in gens(ring)]
        map(ps -> set_prec!(ps, 2 * cur_prec), eval_point)
        eqs_eval = map(p -> evaluate(p, eval_point), eqs)
        J_eval = map(p -> evaluate(p, eval_point), Jac_xs)
        Jd_eval = map(p -> evaluate(p, eval_point), Jac_dots)
        Jd_inv = ps_matrix_inv(Jd_eval, cur_prec)

        X_err = ps_matrix_linear_de(-Jd_inv * J_eval, -Jd_inv * eqs_eval, zero(Svconst))
        for i in 1:length(x_vars)
            solution[x_vars[i]] = solution[x_vars[i]] + X_err[i]
            solution[x_dot_vars[i]] = ps_diff(solution[x_vars[i]])
        end
        cur_prec *= 2
    end

    return solution
end

#------------------------------------------------------------------------------

function ps_ode_solution(
        equations::Array{P, 1}, 
        ic::Dict{P, Int}, 
        inputs::Dict{P, Array{Int, 1}}, 
        prec::Int
    ) where P <: MPolyElem{<: Generic.FieldElem}
    bring = base_ring(parent(equations[1]))
    ps_ode_solution(
        equations,
        Dict(x => bring(v) for (x, v) in ic),
        Dict(u => map(v -> bring(v), vv) for (u, vv) in inputs),
        prec
    )
end
