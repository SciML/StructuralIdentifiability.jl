function truncate_matrix(M, prec)
    return map(e -> truncate(e, prec), M)
end

function fake_truncate_matrix(M, prec)
    prec = precision(M[1, 1])
    truncation = truncate_matrix(M, prec)
    map(e -> set_prec!(e, prec), truncation)
    return truncation
end

function ps_matrix_const_term(M)
    return map(e -> coeff(e, 0), M)
end

function ps_matrix_inv(M)
    """
    Input:
        - M - a square matrix with entries in a univariate power series ring
          it is assumed that M(0) is invertible and all entries having the same precision
    Output: the inverse of M computed up to the precision of the matrix
    """
    const_term = ps_matrix_const_term(M)
    prec = precision(M[1, 1])
    power_series_ring = base_ring(parent(M))
    result = map(a -> power_series_ring(a), Base.inv(const_term))
    cur_prec = 1
    while cur_prec < prec
        result = 2 * result - result * M * result
        cur_prec *= 2
    end
    return result
end

function ps_diff(ps)
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

function ps_integrate(ps)
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

function ps_matrix_log(M)
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

function ps_matrix_homlinear_de(A, Y0)
    """
    Input:
        - A - a square matrix with entries in a univariate power series ring
        - Y0 - a square invertible matrix over the base field
    Output: matrix Y such that Y' = AY up to precision of A - 1 and Y(0) = Y0
    """
    prec = precision(A[1, 1])
    identity = one(parent(A))
    ps_ring = base_ring(parent(A))
    Y = (identity + gen(ps_ring) * fake_truncate_matrix(A, 1)) * map(e -> truncate(ps_ring(e), prec), Y0)
    Z = map(e -> truncate(ps_ring(e), prec), Base.inv(Y0))
    cur_prec = 1
    while cur_prec <= (prec + 1) // 2
        Yprime = map(ps_diff, Y)
        Z = Z + fake_truncate_matrix(Z * (identity - Y * Z), cur_prec)
        Y = Y - fake_truncate_matrix(Y * map(ps_integrate, Z * (Yprime - fake_truncate_matrix(A, 2 * cur_prec - 1) * Y)), 2 * cur_prec)
        cur_prec *= 2
    end
    return Y, Z
end

function ps_matrix_linear_de(A, B, Y0)
    """
    Input:
        - A, B - square matrices with entries in a univariate power series ring
        - Y0 - a square invertible matrix over the base field
    Output: matrix Y such that Y' = AY + B up to precision of A - 1 and Y(0) = Y0
    """
    prec = precision(A[1, 1])
    Yh, Zh = ps_matrix_homlinear_de(A, Y0)
    Zh += fake_truncate_matrix(Zh * (one(parent(A)) - Yh * Zh), prec)
    Y = fake_truncate_matrix(Yh * map(ps_integrate, Zh * B), prec)
    return Y + Yh
end
